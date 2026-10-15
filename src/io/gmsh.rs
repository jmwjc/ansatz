use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::path::Path;

use crate::element::poi1::Poi1;
use crate::element::quad4::Quad4;
use crate::element::seg2::Seg2;
use crate::element::tri3::Tri3;
use crate::node::Node;

// ============================================================================
// Public API
// ============================================================================

/// Bridge trait between an element type and the GMSH file format.
///
/// Element types implementing this trait can be constructed from the
/// intermediate representation as a `Vec<T>` by physical name via
/// [`GmshMesh::elements`]. To support a new element type, just implement
/// this trait for it in the io layer.
pub trait FromGmsh<const D: usize>: Sized {
    /// GMSH v4 element type number (the 3rd field of each block in `$Elements`).
    const GMSH_ELEMENT_TYPE: i32;

    /// Constructs an element from a list of node ids; returns `None` if any node is missing.
    fn from_gmsh(node_ids: &[usize], nodes: &HashMap<usize, Node<D>>) -> Option<Self>;
}

impl<const D: usize> FromGmsh<D> for Seg2<D> {
    const GMSH_ELEMENT_TYPE: i32 = 1; // 2-node line

    fn from_gmsh(node_ids: &[usize], nodes: &HashMap<usize, Node<D>>) -> Option<Self> {
        let n1 = nodes.get(node_ids.first()?)?;
        let n2 = nodes.get(node_ids.get(1)?)?;
        Some(Seg2::new(n1.clone(), n2.clone()))
    }
}

impl<const D: usize> FromGmsh<D> for Tri3<D> {
    const GMSH_ELEMENT_TYPE: i32 = 2; // 3-node triangle

    fn from_gmsh(node_ids: &[usize], nodes: &HashMap<usize, Node<D>>) -> Option<Self> {
        let n1 = nodes.get(node_ids.first()?)?;
        let n2 = nodes.get(node_ids.get(1)?)?;
        let n3 = nodes.get(node_ids.get(2)?)?;
        Some(Tri3::new(n1.clone(), n2.clone(), n3.clone()))
    }
}

impl<const D: usize> FromGmsh<D> for Quad4<D> {
    const GMSH_ELEMENT_TYPE: i32 = 3; // 4-node quadrangle

    fn from_gmsh(node_ids: &[usize], nodes: &HashMap<usize, Node<D>>) -> Option<Self> {
        let n1 = nodes.get(node_ids.first()?)?;
        let n2 = nodes.get(node_ids.get(1)?)?;
        let n3 = nodes.get(node_ids.get(2)?)?;
        let n4 = nodes.get(node_ids.get(3)?)?;
        Some(Quad4::new(n1.clone(), n2.clone(), n3.clone(), n4.clone()))
    }
}

impl<const D: usize> FromGmsh<D> for Poi1<D> {
    const GMSH_ELEMENT_TYPE: i32 = 15; // 1-node point

    fn from_gmsh(node_ids: &[usize], nodes: &HashMap<usize, Node<D>>) -> Option<Self> {
        let n1 = nodes.get(node_ids.first()?)?;
        Some(Poi1::new(n1.clone()))
    }
}

/// Element block of a single physical group: element type number + per-element node id lists.
struct ElementBlock {
    elem_type: i32,
    connectivity: Vec<Vec<usize>>,
}

/// Intermediate product of a GMSH mesh.
///
/// After [`read_gmsh`] parses an `.msh` file, the node coordinate table and the
/// element connectivity grouped by physical group are stored here, without
/// binding to a concrete element type. Callers obtain a typed element vector
/// `Vec<T>` for a physical name via [`GmshMesh::elements`]; since all elements
/// in one physical group share the same element type, `T` in `Vec<T>` is a
/// single type.
pub struct GmshMesh<const D: usize> {
    /// Node table: `node-id → Node`.
    pub nodes: HashMap<usize, Node<D>>,
    /// Physical group names: `physical-tag → name`.
    pub physical_names: HashMap<i32, String>,
    /// Element blocks: `physical-tag → ElementBlock`.
    blocks: HashMap<i32, ElementBlock>,
}

impl<const D: usize> GmshMesh<D> {
    /// Builds the typed element vector for a physical name.
    ///
    /// # Panics
    ///
    /// Panics if the physical name is not found, or if the element type of that
    /// physical group does not match `T`.
    pub fn elements<T: FromGmsh<D>>(&self, physical_name: &str) -> Vec<T> {
        let tag = self
            .physical_names
            .iter()
            .find(|(_, name)| name.as_str() == physical_name)
            .map(|(tag, _)| *tag)
            .unwrap_or_else(|| panic!("physical name not found: {physical_name}"));
        self.elements_by_tag::<T>(tag)
    }

    /// Builds the typed element vector for a physical tag.
    ///
    /// # Panics
    ///
    /// Panics if the physical tag is not found, or if the element type of that
    /// physical group does not match `T`.
    pub fn elements_by_tag<T: FromGmsh<D>>(&self, tag: i32) -> Vec<T> {
        let block = self
            .blocks
            .get(&tag)
            .unwrap_or_else(|| panic!("physical group {tag} has no elements"));
        assert_eq!(
            block.elem_type,
            T::GMSH_ELEMENT_TYPE,
            "physical group {tag} holds element type {}, not {}",
            block.elem_type,
            T::GMSH_ELEMENT_TYPE
        );
        block
            .connectivity
            .iter()
            .filter_map(|ids| T::from_gmsh(ids, &self.nodes))
            .collect()
    }
}

/// Parses a GMSH v4 `.msh` file and returns the intermediate product [`GmshMesh`].
///
/// This function only reads the file and does not construct concrete element
/// types; element vectors are materialized on demand by [`GmshMesh::elements`]
/// per physical name.
///
/// # Panics
///
/// Panics if the file is not GMSH v4, is malformed, or lacks a required section.
pub fn read_gmsh<const D: usize>(filepath: impl AsRef<Path>) -> GmshMesh<D> {
    let file = File::open(filepath.as_ref()).expect("cannot open .msh file");
    let mut reader = BufReader::with_capacity(64 * 1024, file); // 64 KiB buffer

    // ---------- $MeshFormat ----------
    skip_to_section(&mut reader, "$MeshFormat");
    let (version, _binary, _data_size) = read_mesh_format(&mut reader);
    assert!(
        (4.0..5.0).contains(&version),
        "only GMSH v4 format is supported (got {version})"
    );

    // ---------- $PhysicalNames (optional) ----------
    let physical_names: HashMap<i32, String> =
        if peek_section(&mut reader).as_deref() == Some("$PhysicalNames") {
            skip_to_section(&mut reader, "$PhysicalNames");
            read_physical_names(&mut reader)
        } else {
            HashMap::new()
        };

    // ---------- $Entities ----------
    let entity_to_phys: HashMap<(i32, i32), i32> =
        if peek_section(&mut reader).as_deref() == Some("$Entities") {
            skip_to_section(&mut reader, "$Entities");
            read_entities(&mut reader)
        } else {
            panic!("$Entities section not found – only GMSH v4 is supported");
        };

    // ---------- $Nodes ----------
    skip_to_section(&mut reader, "$Nodes");
    let nodes = read_nodes_v4::<D>(&mut reader);

    // ---------- $Elements ----------
    skip_to_section(&mut reader, "$Elements");
    let blocks = read_elements_v4(&mut reader, &entity_to_phys);

    GmshMesh {
        nodes,
        physical_names,
        blocks,
    }
}

// ============================================================================
// Internal helpers – section navigation
// ============================================================================

/// Advance the reader past the line that contains `$SectionName` (inclusive).
fn skip_to_section(reader: &mut BufReader<File>, section: &str) {
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line).expect("I/O error reading .msh");
        assert_ne!(n, 0, "unexpected EOF while looking for {section}");
        if line.trim() == section {
            return;
        }
    }
}

/// Peek at the next non-empty section tag without consuming it.
fn peek_section(reader: &mut BufReader<File>) -> Option<String> {
    let pos = reader.stream_position().ok()?;
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line).ok()?;
        if n == 0 {
            break;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        reader.seek(SeekFrom::Start(pos)).ok()?;
        return Some(trimmed.to_string());
    }
    reader.seek(SeekFrom::Start(pos)).ok()?;
    None
}

// ============================================================================
// $MeshFormat
// ============================================================================

fn read_mesh_format(reader: &mut BufReader<File>) -> (f64, bool, usize) {
    let mut line = String::new();
    reader
        .read_line(&mut line)
        .expect("I/O error reading $MeshFormat");
    let parts: Vec<&str> = line.split_whitespace().collect();
    assert!(parts.len() >= 3, "malformed $MeshFormat line");
    let version: f64 = parts[0].parse().expect("invalid mesh format version");
    let binary: bool = parts[1] != "0";
    let _data_size: usize = parts[2].parse().unwrap_or(8);
    assert!(!binary, "binary GMSH files are not supported");
    let mut end_line = String::new();
    reader.read_line(&mut end_line).expect("I/O error");
    assert!(
        end_line.trim() == "$EndMeshFormat",
        "expected $EndMeshFormat"
    );
    (version, binary, _data_size)
}

// ============================================================================
// $PhysicalNames
// ============================================================================

fn read_physical_names(reader: &mut BufReader<File>) -> HashMap<i32, String> {
    let mut line = String::new();
    reader.read_line(&mut line).expect("I/O error");
    let count: usize = line.trim().parse().expect("invalid $PhysicalNames count");
    let mut map = HashMap::with_capacity(count);
    for _ in 0..count {
        line.clear();
        reader.read_line(&mut line).expect("I/O error");
        let parts: Vec<&str> = line.split_whitespace().collect();
        assert!(parts.len() >= 3, "malformed $PhysicalNames entry");
        let phys_tag: i32 = parts[1].parse().expect("invalid physical tag");
        // Re-join everything after the first two tokens to get the name
        let name_start = parts[0].len() + parts[1].len() + 2;
        let name = line[name_start..].trim().trim_matches('"').to_string();
        map.insert(phys_tag, name);
    }
    let mut end_line = String::new();
    reader.read_line(&mut end_line).expect("I/O error");
    assert!(
        end_line.trim() == "$EndPhysicalNames",
        "expected $EndPhysicalNames"
    );
    map
}

// ============================================================================
// $Entities
// ============================================================================

/// Build `(dim, entity_tag) → physical_group` mapping.
fn read_entities(reader: &mut BufReader<File>) -> HashMap<(i32, i32), i32> {
    let mut line = String::new();
    reader.read_line(&mut line).expect("I/O error");
    let counts: Vec<usize> = line
        .split_whitespace()
        .map(|s| s.parse().expect("invalid $Entities counts"))
        .collect();
    assert_eq!(counts.len(), 4, "expected 4 entity counts");

    let total = counts.iter().sum::<usize>();
    let mut map = HashMap::with_capacity(total);

    read_entity_block(reader, &mut map, 0, counts[0], 3);
    read_entity_block(reader, &mut map, 1, counts[1], 6);
    read_entity_block(reader, &mut map, 2, counts[2], 6);
    read_entity_block(reader, &mut map, 3, counts[3], 6);

    let mut end_line = String::new();
    reader.read_line(&mut end_line).expect("I/O error");
    assert!(end_line.trim() == "$EndEntities", "expected $EndEntities");

    map
}

fn read_entity_block(
    reader: &mut BufReader<File>,
    map: &mut HashMap<(i32, i32), i32>,
    dim: i32,
    count: usize,
    coord_fields: usize,
) {
    let mut line = String::new();
    for _ in 0..count {
        line.clear();
        reader.read_line(&mut line).expect("I/O error");
        let tokens: Vec<&str> = line.split_whitespace().collect();
        let tag: i32 = tokens[0].parse().expect("invalid entity tag");
        let idx = 1 + coord_fields;
        let num_phys: usize = tokens[idx].parse().expect("invalid numPhysicalTags");
        let phys_tag: i32 = if num_phys > 0 {
            tokens[idx + 1].parse().expect("invalid physicalTag")
        } else {
            0
        };
        map.insert((dim, tag), phys_tag);
    }
}

// ============================================================================
// $Nodes (v4)
// ============================================================================

fn read_nodes_v4<const D: usize>(
    reader: &mut BufReader<File>,
) -> HashMap<usize, crate::node::Node<D>> {
    let mut line = String::new();
    reader.read_line(&mut line).expect("I/O error");
    let header: Vec<usize> = line
        .split_whitespace()
        .map(|s| s.parse().expect("invalid $Nodes header"))
        .collect();
    assert_eq!(header.len(), 4, "expected 4 $Nodes header fields");
    let num_entity_blocks = header[0];
    let num_nodes = header[1];

    let mut nodes = HashMap::with_capacity(num_nodes);

    for _ in 0..num_entity_blocks {
        line.clear();
        reader.read_line(&mut line).expect("I/O error");
        let block_header: Vec<&str> = line.split_whitespace().collect();
        assert!(
            block_header.len() >= 4,
            "malformed node entity block header"
        );
        let _entity_dim: i32 = block_header[0].parse().expect("invalid entityDim");
        let parametric: i32 = block_header[2].parse().expect("invalid parametric flag");
        let num_nodes_in_block: usize = block_header[3].parse().expect("invalid node count");

        // Auto-detect: read the first data line to determine format.
        line.clear();
        reader.read_line(&mut line).expect("I/O error");
        let first_tokens: Vec<&str> = line.split_whitespace().collect();

        // If first line has >1 tokens, it's compact format: "tag X Y Z" per line.
        // Otherwise (1 token) it's separated: tags first, then coordinates.
        let separated = first_tokens.len() == 1;

        assert_eq!(
            parametric, 0,
            "parametric nodes only supported in separated format"
        );

        if separated {
            // Separated format: tags then coordinates.
            // Re-process the first line which contains the first tag(s).
            let mut tags: Vec<usize> = Vec::with_capacity(num_nodes_in_block);
            let first_tags: Vec<&str> = first_tokens;
            for t in &first_tags {
                if tags.len() >= num_nodes_in_block {
                    break;
                }
                if !t.is_empty() {
                    tags.push(t.parse().expect("invalid node tag"));
                }
            }
            while tags.len() < num_nodes_in_block {
                line.clear();
                reader.read_line(&mut line).expect("I/O error");
                let toks: Vec<&str> = line.split_whitespace().collect();
                for t in &toks {
                    if tags.len() >= num_nodes_in_block {
                        break;
                    }
                    if !t.is_empty() {
                        tags.push(t.parse().expect("invalid node tag"));
                    }
                }
            }
            // Now read coordinates.
            for tag in tags.iter().take(num_nodes_in_block) {
                line.clear();
                reader.read_line(&mut line).expect("I/O error");
                let toks: Vec<&str> = line.split_whitespace().collect();
                let ncoords = D.min(toks.len());
                let mut coords = [0.0_f64; D];
                for d in 0..ncoords {
                    coords[d] = toks[d].parse().expect("invalid node coordinate");
                }
                nodes.insert(*tag - 1, crate::node::Node::new(*tag - 1, coords));
            }
        } else {
            // Compact format: each line is "tag x y z".
            let first_id: usize = first_tokens[0].parse().expect("invalid node tag");
            let ncoords = D.min(first_tokens.len() - 1);
            let mut first_coords = [0.0_f64; D];
            for d in 0..ncoords {
                first_coords[d] = first_tokens[1 + d]
                    .parse()
                    .expect("invalid node coordinate");
            }
            nodes.insert(first_id - 1, crate::node::Node::new(first_id - 1, first_coords));

            for _node_idx in 1..num_nodes_in_block {
                line.clear();
                reader.read_line(&mut line).expect("I/O error");
                let toks: Vec<&str> = line.split_whitespace().collect();
                let id: usize = toks[0].parse().expect("invalid node tag");
                let n = D.min(toks.len() - 1);
                let mut coords = [0.0_f64; D];
                for d in 0..n {
                    coords[d] = toks[1 + d].parse().expect("invalid node coordinate");
                }
                nodes.insert(id - 1, crate::node::Node::new(id - 1, coords));
            }
        }
    }

    let mut end_line = String::new();
    reader.read_line(&mut end_line).expect("I/O error");
    assert!(end_line.trim() == "$EndNodes", "expected $EndNodes");

    nodes
}

// ============================================================================
// $Elements (v4)
// ============================================================================

/// Read raw element connectivity grouped by physical group, without constructing
/// concrete element types. Each block's element type and per-element node ids are
/// preserved so callers can materialize typed vectors on demand.
fn read_elements_v4(
    reader: &mut BufReader<File>,
    entity_to_phys: &HashMap<(i32, i32), i32>,
) -> HashMap<i32, ElementBlock> {
    let mut line = String::new();
    reader.read_line(&mut line).expect("I/O error");
    let header: Vec<usize> = line
        .split_whitespace()
        .map(|s| s.parse().expect("invalid $Elements header"))
        .collect();
    assert_eq!(header.len(), 4, "expected 4 $Elements header fields");
    let num_entity_blocks = header[0];

    let mut blocks: HashMap<i32, ElementBlock> = HashMap::new();

    for _block in 0..num_entity_blocks {
        line.clear();
        reader.read_line(&mut line).expect("I/O error");
        let block_header: Vec<&str> = line.split_whitespace().collect();
        assert!(
            block_header.len() >= 4,
            "malformed element entity block header"
        );
        let entity_dim: i32 = block_header[0].parse().expect("invalid entityDim");
        let entity_tag: i32 = block_header[1].parse().expect("invalid entityTag");
        let elem_type: i32 = block_header[2].parse().expect("invalid elementType");
        let num_elems_in_block: usize = block_header[3].parse().expect("invalid element count");

        let phys_group = entity_to_phys
            .get(&(entity_dim, entity_tag))
            .copied()
            .unwrap_or(0);

        let block = blocks
            .entry(phys_group)
            .or_insert_with(|| ElementBlock {
                elem_type,
                connectivity: Vec::with_capacity(num_elems_in_block),
            });
        assert_eq!(
            block.elem_type, elem_type,
            "physical group {phys_group} mixes element types {} and {elem_type}",
            block.elem_type
        );
        block.connectivity.reserve(num_elems_in_block);

        for _elem in 0..num_elems_in_block {
            line.clear();
            reader.read_line(&mut line).expect("I/O error");
            let tokens: Vec<&str> = line.split_whitespace().collect();
            // elementTag nodeTag1 nodeTag2 ...
            let node_ids: Vec<usize> = tokens[1..]
                .iter()
                .map(|s| s.parse::<usize>().expect("invalid node tag") - 1)
                .collect();
            block.connectivity.push(node_ids);
        }
    }

    let mut end_line = String::new();
    reader.read_line(&mut end_line).expect("I/O error");
    assert!(end_line.trim() == "$EndElements", "expected $EndElements");

    for b in blocks.values_mut() {
        b.connectivity.shrink_to_fit();
    }

    blocks
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::Element;

    #[test]
    fn reads_segments_and_triangles() {
        let mesh = read_gmsh::<2>("examples/msh/patchtest_tri3.msh");
        let segs: Vec<Seg2<2>> = mesh.elements("Γᵍ");
        assert_eq!(segs.len(), 40, "expected 40 segment elements");
        let tris: Vec<Tri3<2>> = mesh.elements("Ω");
        assert_eq!(tris.len(), 288, "expected 288 triangle elements");
    }

    #[test]
    fn quad4_from_gmsh() {
        let mut nodes = HashMap::new();
        nodes.insert(0, Node::new(0, [0.0, 0.0]));
        nodes.insert(1, Node::new(1, [1.0, 0.0]));
        nodes.insert(2, Node::new(2, [1.0, 1.0]));
        nodes.insert(3, Node::new(3, [0.0, 1.0]));
        let q = Quad4::from_gmsh(&[0, 1, 2, 3], &nodes).unwrap();
        assert_eq!(q.id(), [0, 1, 2, 3]);
    }

    #[test]
    fn reads_points_and_segments() {
        let mesh = read_gmsh::<1>("examples/msh/bar1d.msh");
        // phys 0 (unnamed) → 2 point elements; phys 1 "line" → 1 segment element.
        let points: Vec<Poi1<1>> = mesh.elements_by_tag(0);
        assert_eq!(points.len(), 2, "expected 2 point elements");
        let segs: Vec<Seg2<1>> = mesh.elements("line");
        assert_eq!(segs.len(), 1, "expected 1 segment element");
    }
}