use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use crate::element::poi1::Poi1;
use crate::element::quad4::Quad4;
use crate::element::seg2::Seg2;
use crate::element::tri3::Tri3;
use crate::node::Node;

// ============================================================================
// Errors
// ============================================================================

/// Errors produced while reading or materializing a GMSH mesh.
#[derive(Debug)]
pub enum MeshError {
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    UnsupportedVersion {
        found: f64,
    },
    BinaryNotSupported,
    MissingSection {
        section: &'static str,
    },
    UnexpectedEof {
        section: &'static str,
    },
    ExpectedSection {
        expected: &'static str,
        found: String,
    },
    MalformedLine {
        section: &'static str,
        line: usize,
        message: String,
    },
    Parse {
        section: &'static str,
        line: usize,
        field: &'static str,
        value: String,
    },
    InvalidNodeTag {
        tag: usize,
        line: usize,
    },
    DuplicateNode {
        tag: usize,
        line: usize,
    },
    NodeCountMismatch {
        expected: usize,
        found: usize,
    },
    ElementCountMismatch {
        expected: usize,
        found: usize,
    },
    UnsupportedParametricNodes {
        line: usize,
    },
    UnknownPhysicalName {
        name: String,
    },
    EmptyPhysicalGroup {
        tag: i32,
    },
    ElementTypeMismatch {
        group: i32,
        expected: i32,
        found: i32,
    },
    InvalidConnectivity {
        group: i32,
        element: usize,
        expected: usize,
        found: usize,
    },
    MissingNode {
        group: i32,
        element: usize,
        tag: usize,
    },
    MissingNodeTag {
        tag: usize,
    },
    MalformedConnectivity {
        expected: usize,
        found: usize,
    },
    NodeIndexOutOfRange {
        index: usize,
        len: usize,
    },
    NodeIdMismatch {
        index: usize,
        id: usize,
    },
    UnknownPhysicalTag {
        tag: i32,
    },
    DuplicatePhysicalName {
        name: String,
    },
}

impl fmt::Display for MeshError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io { path, source } => {
                write!(f, "failed to read {}: {source}", path.display())
            }
            Self::UnsupportedVersion { found } => {
                write!(
                    f,
                    "unsupported GMSH version {found}; only v4 ASCII is supported"
                )
            }
            Self::BinaryNotSupported => {
                write!(f, "binary GMSH files are not supported; export ASCII")
            }
            Self::MissingSection { section } => {
                write!(f, "missing GMSH section {section}")
            }
            Self::UnexpectedEof { section } => {
                write!(f, "unexpected EOF in {section}")
            }
            Self::ExpectedSection { expected, found } => {
                write!(f, "expected {expected}, found {found}")
            }
            Self::MalformedLine {
                section,
                line,
                message,
            } => {
                write!(f, "malformed line {line} in {section}: {message}")
            }
            Self::Parse {
                section,
                line,
                field,
                value,
            } => {
                write!(f, "cannot parse {field} '{value}' at {section}:{line}")
            }
            Self::InvalidNodeTag { tag, line } => {
                write!(
                    f,
                    "invalid node tag {tag} at line {line}; GMSH tags are 1-based"
                )
            }
            Self::DuplicateNode { tag, line } => {
                write!(f, "duplicate node tag {tag} at line {line}")
            }
            Self::NodeCountMismatch { expected, found } => {
                write!(f, "expected {expected} nodes, read {found}")
            }
            Self::ElementCountMismatch { expected, found } => {
                write!(f, "expected {expected} elements, read {found}")
            }
            Self::UnsupportedParametricNodes { line } => {
                write!(f, "parametric nodes are not supported (line {line})")
            }
            Self::UnknownPhysicalName { name } => {
                write!(f, "physical name not found: {name}")
            }
            Self::EmptyPhysicalGroup { tag } => {
                write!(f, "physical group {tag} has no elements")
            }
            Self::ElementTypeMismatch {
                group,
                expected,
                found,
            } => {
                write!(
                    f,
                    "physical group {group} holds element type {found}, expected {expected}"
                )
            }
            Self::InvalidConnectivity {
                group,
                element,
                expected,
                found,
            } => {
                write!(
                    f,
                    "physical group {group}, element {element}: expected {expected} nodes, found {found}"
                )
            }
            Self::MissingNode {
                group,
                element,
                tag,
            } => {
                write!(
                    f,
                    "physical group {group}, element {element} references missing node tag {tag}"
                )
            }
            Self::MissingNodeTag { tag } => {
                write!(f, "node tag {tag} is not present in the node table")
            }
            Self::MalformedConnectivity { expected, found } => {
                write!(
                    f,
                    "malformed connectivity: expected at least {expected} nodes, found {found}"
                )
            }
            Self::NodeIndexOutOfRange { index, len } => {
                write!(
                    f,
                    "node index {index} out of range for node table of length {len}"
                )
            }
            Self::NodeIdMismatch { index, id } => {
                write!(f, "node at index {index} has id {id}; expected {index}")
            }
            Self::UnknownPhysicalTag { tag } => {
                write!(f, "physical tag {tag} is not present")
            }
            Self::DuplicatePhysicalName { name } => {
                write!(f, "duplicate physical name: {name}")
            }
        }
    }
}

impl std::error::Error for MeshError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

// ============================================================================
// Public API
// ============================================================================

/// Dense node coordinate table: node id is the index into `coords`.
///
/// The table is built after the preprocessing renumbering step, so it contains
/// no gaps and no redundant `Node::id` field.
#[derive(Clone, Debug, Default)]
pub struct NodeTable<const D: usize> {
    coords: Vec<[f64; D]>,
}

impl<const D: usize> NodeTable<D> {
    /// Builds a table directly from dense coordinates (`id == index`).
    pub fn new(coords: Vec<[f64; D]>) -> Self {
        Self { coords }
    }

    /// Builds a table from a vector of nodes whose ids are `0..nodes.len()`.
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::NodeIdMismatch`] if a node id does not equal its
    /// position in the vector.
    pub fn from_nodes(nodes: Vec<Node<D>>) -> Result<Self, MeshError> {
        let mut coords = Vec::with_capacity(nodes.len());
        for (index, node) in nodes.into_iter().enumerate() {
            if node.id != index {
                return Err(MeshError::NodeIdMismatch { index, id: node.id });
            }
            coords.push(node.coordinates);
        }
        Ok(Self { coords })
    }

    /// Builds a dense table from a raw node map.
    ///
    /// This is the compatibility path for [`GmshMesh`], whose parser stores
    /// nodes by their raw (tag - 1) ids. The ids must be exactly `0..len`.
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::NodeIndexOutOfRange`] for a non-contiguous id, or
    /// [`MeshError::MissingNodeTag`] if an id in the range is absent.
    pub fn from_map(nodes: &HashMap<usize, Node<D>>) -> Result<Self, MeshError> {
        let len = nodes.len();
        let mut coords = vec![[0.0_f64; D]; len];
        let mut seen = vec![false; len];

        for (&id, node) in nodes {
            if id >= len {
                return Err(MeshError::NodeIndexOutOfRange { index: id, len });
            }
            coords[id] = node.coordinates;
            seen[id] = true;
        }

        if let Some(index) = seen.iter().position(|seen| !seen) {
            return Err(MeshError::MissingNodeTag { tag: index });
        }

        Ok(Self { coords })
    }

    /// Number of nodes.
    pub fn len(&self) -> usize {
        self.coords.len()
    }

    /// Whether the table is empty.
    pub fn is_empty(&self) -> bool {
        self.coords.is_empty()
    }

    /// Returns the coordinates of node `id`.
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::NodeIndexOutOfRange`] if `id` is out of bounds.
    pub fn get(&self, id: usize) -> Result<[f64; D], MeshError> {
        self.coords
            .get(id)
            .copied()
            .ok_or(MeshError::NodeIndexOutOfRange {
                index: id,
                len: self.coords.len(),
            })
    }

    /// All node coordinates, indexed by dense node id.
    pub fn coords(&self) -> &[[f64; D]] {
        &self.coords
    }
}

/// Read-only view of one physical group.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PhysicalGroup<'a> {
    /// GMSH physical tag.
    pub tag: i32,
    /// Physical name, or `None` for an unnamed group.
    pub name: Option<&'a str>,
    /// GMSH element type number; `0` when the group has no element block.
    pub element_type: i32,
    /// Number of elements in the group.
    pub element_count: usize,
}

/// Bridge trait between an element type and the GMSH file format.
///
/// Element types implementing this trait can be constructed from the
/// intermediate representation as a `Vec<T>` by physical name via
/// [`GmshMesh::elements`]. To support a new element type, just implement
/// this trait for it in the io layer.
pub trait FromGmsh<const D: usize>: Sized {
    /// GMSH v4 element type number (the 3rd field of each block in `$Elements`).
    const GMSH_ELEMENT_TYPE: i32;

    /// Number of nodes referenced by one element of this type.
    const GMSH_NODE_COUNT: usize;

    /// Constructs an element from a list of node ids and a dense node table.
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::MissingNodeTag`] if a referenced node is absent,
    /// [`MeshError::NodeIndexOutOfRange`] if an id is outside the table, and
    /// [`MeshError::MalformedConnectivity`] if `node_ids` is too short.
    fn from_gmsh(node_ids: &[usize], nodes: &NodeTable<D>) -> Result<Self, MeshError>;
}

fn node_at<const D: usize>(
    node_ids: &[usize],
    index: usize,
    nodes: &NodeTable<D>,
) -> Result<Node<D>, MeshError> {
    let id = *node_ids
        .get(index)
        .ok_or(MeshError::MalformedConnectivity {
            expected: index + 1,
            found: node_ids.len(),
        })?;
    let coords = nodes.get(id)?;
    Ok(Node::new(id, coords))
}

impl<const D: usize> FromGmsh<D> for Seg2<D> {
    const GMSH_ELEMENT_TYPE: i32 = 1; // 2-node line
    const GMSH_NODE_COUNT: usize = 2;

    fn from_gmsh(node_ids: &[usize], nodes: &NodeTable<D>) -> Result<Self, MeshError> {
        let n1 = node_at(node_ids, 0, nodes)?;
        let n2 = node_at(node_ids, 1, nodes)?;
        Ok(Seg2::new(n1, n2))
    }
}

impl<const D: usize> FromGmsh<D> for Tri3<D> {
    const GMSH_ELEMENT_TYPE: i32 = 2; // 3-node triangle
    const GMSH_NODE_COUNT: usize = 3;

    fn from_gmsh(node_ids: &[usize], nodes: &NodeTable<D>) -> Result<Self, MeshError> {
        let n1 = node_at(node_ids, 0, nodes)?;
        let n2 = node_at(node_ids, 1, nodes)?;
        let n3 = node_at(node_ids, 2, nodes)?;
        Ok(Tri3::new(n1, n2, n3))
    }
}

impl<const D: usize> FromGmsh<D> for Quad4<D> {
    const GMSH_ELEMENT_TYPE: i32 = 3; // 4-node quadrangle
    const GMSH_NODE_COUNT: usize = 4;

    fn from_gmsh(node_ids: &[usize], nodes: &NodeTable<D>) -> Result<Self, MeshError> {
        let n1 = node_at(node_ids, 0, nodes)?;
        let n2 = node_at(node_ids, 1, nodes)?;
        let n3 = node_at(node_ids, 2, nodes)?;
        let n4 = node_at(node_ids, 3, nodes)?;
        Ok(Quad4::new(n1, n2, n3, n4))
    }
}

impl<const D: usize> FromGmsh<D> for Poi1<D> {
    const GMSH_ELEMENT_TYPE: i32 = 15; // 1-node point
    const GMSH_NODE_COUNT: usize = 1;

    fn from_gmsh(node_ids: &[usize], nodes: &NodeTable<D>) -> Result<Self, MeshError> {
        let n1 = node_at(node_ids, 0, nodes)?;
        Ok(Poi1::new(n1))
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
    /// Reverse lookup for names: `name → physical-tag`.
    name_to_tag: HashMap<String, i32>,
    /// Element blocks: `physical-tag → ElementBlock`.
    blocks: HashMap<i32, ElementBlock>,
}

impl<const D: usize> fmt::Debug for GmshMesh<D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut groups: Vec<(i32, i32, usize)> = self
            .blocks
            .iter()
            .map(|(&tag, block)| (tag, block.elem_type, block.connectivity.len()))
            .collect();
        groups.sort_unstable_by_key(|(tag, _, _)| *tag);

        f.debug_struct("GmshMesh")
            .field("n_nodes", &self.nodes.len())
            .field("physical_names", &self.physical_names)
            .field("groups", &groups)
            .finish()
    }
}

impl<const D: usize> GmshMesh<D> {
    /// Resolves a physical name to its physical tag.
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::UnknownPhysicalName`] if the name is not present.
    fn physical_tag(&self, physical_name: &str) -> Result<i32, MeshError> {
        if let Some(&tag) = self.name_to_tag.get(physical_name) {
            return Ok(tag);
        }

        // Fallback for callers that mutate the public `physical_names` map
        // after construction.
        self.physical_names
            .iter()
            .find(|(_, name)| name.as_str() == physical_name)
            .map(|(tag, _)| *tag)
            .ok_or_else(|| MeshError::UnknownPhysicalName {
                name: physical_name.to_string(),
            })
    }

    /// Returns all physical groups in ascending tag order.
    ///
    /// The iterator includes named groups without elements (`element_type == 0`,
    /// `element_count == 0`) and unnamed groups that have elements
    /// (`name == None`).
    pub fn physical_groups(&self) -> impl Iterator<Item = PhysicalGroup<'_>> + '_ {
        let mut tags: Vec<i32> = self
            .blocks
            .keys()
            .chain(self.physical_names.keys())
            .copied()
            .collect();
        tags.sort_unstable();
        tags.dedup();

        tags.into_iter()
            .map(move |tag| self.physical_group_data(tag))
    }

    /// Returns one physical group by tag.
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::UnknownPhysicalTag`] if the tag is not present in
    /// either the physical-name table or the element blocks.
    pub fn physical_group(&self, tag: i32) -> Result<PhysicalGroup<'_>, MeshError> {
        if !self.blocks.contains_key(&tag) && !self.physical_names.contains_key(&tag) {
            return Err(MeshError::UnknownPhysicalTag { tag });
        }
        Ok(self.physical_group_data(tag))
    }

    /// Returns one physical group by name.
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::UnknownPhysicalName`] if the name is not present.
    pub fn physical_group_by_name(&self, name: &str) -> Result<PhysicalGroup<'_>, MeshError> {
        let tag = self.physical_tag(name)?;
        Ok(self.physical_group_data(tag))
    }

    /// Returns named physical groups `(tag, name)` in ascending tag order.
    pub fn physical_names(&self) -> impl Iterator<Item = (i32, &str)> + '_ {
        let mut tags: Vec<i32> = self.physical_names.keys().copied().collect();
        tags.sort_unstable();
        tags.into_iter().filter_map(move |tag| {
            self.physical_names
                .get(&tag)
                .map(|name| (tag, name.as_str()))
        })
    }

    /// Whether a physical group with this name exists.
    pub fn has_physical_group(&self, name: &str) -> bool {
        self.physical_tag(name).is_ok()
    }

    fn physical_group_data(&self, tag: i32) -> PhysicalGroup<'_> {
        let block = self.blocks.get(&tag);
        PhysicalGroup {
            tag,
            name: self.physical_names.get(&tag).map(String::as_str),
            element_type: block.map_or(0, |block| block.elem_type),
            element_count: block.map_or(0, |block| block.connectivity.len()),
        }
    }

    /// Builds a dense [`NodeTable`] from the raw node map.
    ///
    /// # Errors
    ///
    /// Returns an error if node ids are not exactly `0..nodes.len()`.
    pub fn node_table(&self) -> Result<NodeTable<D>, MeshError> {
        NodeTable::from_map(&self.nodes)
    }

    /// Returns the GMSH element type number stored in a physical group.
    ///
    /// This lets callers choose a concrete element type (e.g. `Tri3` vs
    /// `Quad4`) at run time without attempting a typed [`GmshMesh::elements`]
    /// call first.
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::UnknownPhysicalName`] if the name is not present,
    /// or [`MeshError::EmptyPhysicalGroup`] if that group has no elements.
    pub fn element_type(&self, physical_name: &str) -> Result<i32, MeshError> {
        let tag = self.physical_tag(physical_name)?;
        self.blocks
            .get(&tag)
            .map(|block| block.elem_type)
            .ok_or(MeshError::EmptyPhysicalGroup { tag })
    }

    /// Builds the typed element vector for a physical name.
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::UnknownPhysicalName`] if the name is not present,
    /// [`MeshError::ElementTypeMismatch`] if the group holds a different
    /// element type, or a connectivity error if an element cannot be built.
    pub fn elements<T: FromGmsh<D>>(&self, physical_name: &str) -> Result<Vec<T>, MeshError> {
        let tag = self.physical_tag(physical_name)?;
        self.elements_by_tag(tag)
    }

    /// Builds the typed element vector for a physical tag.
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::EmptyPhysicalGroup`] if the tag is not present,
    /// [`MeshError::ElementTypeMismatch`] if the group holds a different
    /// element type, or a connectivity error if an element cannot be built.
    pub fn elements_by_tag<T: FromGmsh<D>>(&self, tag: i32) -> Result<Vec<T>, MeshError> {
        let node_table = self.node_table()?;
        self.elements_by_tag_with_table(tag, &node_table)
    }

    /// Builds the typed element vector for a physical tag using a caller-provided
    /// dense [`NodeTable`].
    ///
    /// This is the reuse-friendly variant of [`GmshMesh::elements_by_tag`]:
    /// build the table once with [`GmshMesh::node_table`] and pass it to every
    /// group.
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::EmptyPhysicalGroup`] if the tag is not present,
    /// [`MeshError::ElementTypeMismatch`] if the group holds a different
    /// element type, or a connectivity error if an element cannot be built.
    pub fn elements_by_tag_with_table<T: FromGmsh<D>>(
        &self,
        tag: i32,
        node_table: &NodeTable<D>,
    ) -> Result<Vec<T>, MeshError> {
        let block = self
            .blocks
            .get(&tag)
            .ok_or(MeshError::EmptyPhysicalGroup { tag })?;

        if block.elem_type != T::GMSH_ELEMENT_TYPE {
            return Err(MeshError::ElementTypeMismatch {
                group: tag,
                expected: T::GMSH_ELEMENT_TYPE,
                found: block.elem_type,
            });
        }

        let mut out = Vec::with_capacity(block.connectivity.len());
        for (element, ids) in block.connectivity.iter().enumerate() {
            if ids.len() != T::GMSH_NODE_COUNT {
                return Err(MeshError::InvalidConnectivity {
                    group: tag,
                    element,
                    expected: T::GMSH_NODE_COUNT,
                    found: ids.len(),
                });
            }

            let elm = T::from_gmsh(ids, node_table).map_err(|err| match err {
                MeshError::MissingNodeTag { tag: node_tag } => MeshError::MissingNode {
                    group: tag,
                    element,
                    tag: node_tag,
                },
                MeshError::NodeIndexOutOfRange { index, .. } => MeshError::MissingNode {
                    group: tag,
                    element,
                    tag: index,
                },
                MeshError::MalformedConnectivity { expected, found } => {
                    MeshError::InvalidConnectivity {
                        group: tag,
                        element,
                        expected,
                        found,
                    }
                }
                other => other,
            })?;

            out.push(elm);
        }

        Ok(out)
    }
}

/// Parses a GMSH v4 `.msh` file and returns the intermediate product [`GmshMesh`].
///
/// This function only reads the file and does not construct concrete element
/// types; element vectors are materialized on demand by [`GmshMesh::elements`]
/// per physical name.
///
/// # Errors
///
/// Returns a [`MeshError`] if the file cannot be opened, is not GMSH v4 ASCII,
/// is malformed, or lacks a required section.
pub fn read_gmsh<const D: usize>(filepath: impl AsRef<Path>) -> Result<GmshMesh<D>, MeshError> {
    let mut reader = MshReader::open(filepath.as_ref())?;

    // ---------- $MeshFormat ----------
    reader.expect_section("$MeshFormat")?;
    let _ = read_mesh_format(&mut reader)?;

    // ---------- $PhysicalNames (optional) ----------
    let physical_names: HashMap<i32, String> =
        if reader.peek_section()?.as_deref() == Some("$PhysicalNames") {
            reader.expect_section("$PhysicalNames")?;
            read_physical_names(&mut reader)?
        } else {
            HashMap::new()
        };

    let name_to_tag = build_name_to_tag(&physical_names)?;

    // ---------- $Entities ----------
    let entity_to_phys: HashMap<(i32, i32), i32> =
        if reader.peek_section()?.as_deref() == Some("$Entities") {
            reader.expect_section("$Entities")?;
            read_entities(&mut reader)?
        } else {
            return Err(MeshError::MissingSection {
                section: "$Entities",
            });
        };

    // ---------- $Nodes ----------
    reader.expect_section("$Nodes")?;
    let nodes = read_nodes_v4::<D>(&mut reader)?;

    // ---------- $Elements ----------
    reader.expect_section("$Elements")?;
    let blocks = read_elements_v4(&mut reader, &entity_to_phys)?;

    Ok(GmshMesh {
        nodes,
        physical_names,
        name_to_tag,
        blocks,
    })
}

// ============================================================================
// Reader with line and section context
// ============================================================================

struct MshReader {
    inner: BufReader<File>,
    path: PathBuf,
    line: usize,
    section: &'static str,
}

impl MshReader {
    fn open(path: &Path) -> Result<Self, MeshError> {
        let file = File::open(path).map_err(|source| MeshError::Io {
            path: path.to_owned(),
            source,
        })?;

        Ok(Self {
            inner: BufReader::with_capacity(64 * 1024, file),
            path: path.to_owned(),
            line: 0,
            section: "",
        })
    }

    fn io_error(&self, source: std::io::Error) -> MeshError {
        MeshError::Io {
            path: self.path.clone(),
            source,
        }
    }

    fn read_line(&mut self) -> Result<String, MeshError> {
        let mut line = String::new();
        let n = self
            .inner
            .read_line(&mut line)
            .map_err(|e| self.io_error(e))?;
        if n == 0 {
            return Err(MeshError::UnexpectedEof {
                section: self.section,
            });
        }
        self.line += 1;
        Ok(line)
    }

    fn expect_section(&mut self, section: &'static str) -> Result<(), MeshError> {
        self.section = section;

        loop {
            let line = self.read_line()?;
            let trimmed = line.trim();
            if trimmed == section {
                return Ok(());
            }
            if trimmed.starts_with('$') && !trimmed.starts_with("$End") {
                return Err(MeshError::ExpectedSection {
                    expected: section,
                    found: trimmed.to_string(),
                });
            }
        }
    }

    fn expect_end(&mut self, end: &'static str) -> Result<(), MeshError> {
        let line = self.read_line()?;
        let trimmed = line.trim();
        if trimmed != end {
            return Err(MeshError::ExpectedSection {
                expected: end,
                found: trimmed.to_string(),
            });
        }
        Ok(())
    }

    fn peek_section(&mut self) -> Result<Option<String>, MeshError> {
        let pos = self.inner.stream_position().map_err(|e| self.io_error(e))?;

        loop {
            let mut line = String::new();
            let n = self
                .inner
                .read_line(&mut line)
                .map_err(|e| self.io_error(e))?;
            if n == 0 {
                break;
            }

            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }

            self.inner
                .seek(SeekFrom::Start(pos))
                .map_err(|e| self.io_error(e))?;
            return Ok(Some(trimmed.to_string()));
        }

        self.inner
            .seek(SeekFrom::Start(pos))
            .map_err(|e| self.io_error(e))?;
        Ok(None)
    }

    fn parse<T: std::str::FromStr>(
        &self,
        value: &str,
        field: &'static str,
    ) -> Result<T, MeshError> {
        value.parse().map_err(|_| MeshError::Parse {
            section: self.section,
            line: self.line,
            field,
            value: value.to_string(),
        })
    }

    fn malformed(&self, message: impl Into<String>) -> MeshError {
        MeshError::MalformedLine {
            section: self.section,
            line: self.line,
            message: message.into(),
        }
    }
}

// ============================================================================
// $MeshFormat
// ============================================================================

fn read_mesh_format(reader: &mut MshReader) -> Result<(f64, bool, usize), MeshError> {
    let line = reader.read_line()?;
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 3 {
        return Err(reader.malformed("expected: version binary data_size"));
    }

    let version: f64 = reader.parse(parts[0], "version")?;
    let binary = parts[1] != "0";
    let data_size: usize = reader.parse(parts[2], "data_size")?;

    if !(4.0..5.0).contains(&version) {
        return Err(MeshError::UnsupportedVersion { found: version });
    }
    if binary {
        return Err(MeshError::BinaryNotSupported);
    }

    reader.expect_end("$EndMeshFormat")?;
    Ok((version, binary, data_size))
}

// ============================================================================
// $PhysicalNames
// ============================================================================

fn read_physical_names(reader: &mut MshReader) -> Result<HashMap<i32, String>, MeshError> {
    let count_line = reader.read_line()?;
    let count: usize = reader.parse(count_line.trim(), "physical name count")?;
    let mut map = HashMap::with_capacity(count);

    for _ in 0..count {
        let line = reader.read_line()?;
        let line = line.trim();

        let (dim, rest) = line
            .split_once(|c: char| c.is_whitespace())
            .ok_or_else(|| reader.malformed("expected: dimension tag \"name\""))?;
        let _dim: i32 = reader.parse(dim, "dimension")?;

        let rest = rest.trim_start();
        let (tag, name) = rest
            .split_once(|c: char| c.is_whitespace())
            .ok_or_else(|| reader.malformed("expected: tag \"name\""))?;
        let tag: i32 = reader.parse(tag, "physical tag")?;

        let name = name.trim().trim_matches('"').to_string();
        if map.insert(tag, name).is_some() {
            return Err(reader.malformed(format!("duplicate physical tag {tag}")));
        }
    }

    reader.expect_end("$EndPhysicalNames")?;
    Ok(map)
}

fn build_name_to_tag(names: &HashMap<i32, String>) -> Result<HashMap<String, i32>, MeshError> {
    let mut map = HashMap::with_capacity(names.len());
    for (&tag, name) in names {
        if map.insert(name.clone(), tag).is_some() {
            return Err(MeshError::DuplicatePhysicalName { name: name.clone() });
        }
    }
    Ok(map)
}

// ============================================================================
// $Entities
// ============================================================================

/// Build `(dim, entity_tag) → physical_group` mapping.
fn read_entities(reader: &mut MshReader) -> Result<HashMap<(i32, i32), i32>, MeshError> {
    let counts: Vec<usize> = reader
        .read_line()?
        .split_whitespace()
        .map(|s| reader.parse(s, "entity count"))
        .collect::<Result<_, _>>()?;

    if counts.len() != 4 {
        return Err(reader.malformed("expected 4 entity counts"));
    }

    let total = counts.iter().sum::<usize>();
    let mut map = HashMap::with_capacity(total);

    read_entity_block(reader, &mut map, 0, counts[0], 3)?;
    read_entity_block(reader, &mut map, 1, counts[1], 6)?;
    read_entity_block(reader, &mut map, 2, counts[2], 6)?;
    read_entity_block(reader, &mut map, 3, counts[3], 6)?;

    reader.expect_end("$EndEntities")?;
    Ok(map)
}

fn read_entity_block(
    reader: &mut MshReader,
    map: &mut HashMap<(i32, i32), i32>,
    dim: i32,
    count: usize,
    coord_fields: usize,
) -> Result<(), MeshError> {
    for _ in 0..count {
        let line = reader.read_line()?;
        let tokens: Vec<&str> = line.split_whitespace().collect();

        let tag: i32 = reader.parse(
            tokens
                .first()
                .ok_or_else(|| reader.malformed("missing entity tag"))?,
            "entity tag",
        )?;

        let idx = 1 + coord_fields;
        let num_phys: usize = reader.parse(
            tokens
                .get(idx)
                .ok_or_else(|| reader.malformed("missing numPhysicalTags"))?,
            "numPhysicalTags",
        )?;

        let phys_tag: i32 = if num_phys > 0 {
            reader.parse(
                tokens
                    .get(idx + 1)
                    .ok_or_else(|| reader.malformed("missing physicalTag"))?,
                "physicalTag",
            )?
        } else {
            0
        };

        map.insert((dim, tag), phys_tag);
    }
    Ok(())
}

// ============================================================================
// $Nodes (v4)
// ============================================================================

fn read_nodes_v4<const D: usize>(
    reader: &mut MshReader,
) -> Result<HashMap<usize, Node<D>>, MeshError> {
    let header: Vec<usize> = reader
        .read_line()?
        .split_whitespace()
        .map(|s| reader.parse(s, "node header"))
        .collect::<Result<_, _>>()?;

    if header.len() != 4 {
        return Err(reader.malformed("expected 4 $Nodes header fields"));
    }
    let num_entity_blocks = header[0];
    let num_nodes = header[1];

    let mut nodes = HashMap::with_capacity(num_nodes);

    for _ in 0..num_entity_blocks {
        let block_line = reader.read_line()?;
        let block_header: Vec<&str> = block_line.split_whitespace().collect();
        if block_header.len() < 4 {
            return Err(reader.malformed("malformed node entity block header"));
        }

        let _entity_dim: i32 = reader.parse(block_header[0], "entityDim")?;
        let _entity_tag: i32 = reader.parse(block_header[1], "entityTag")?;
        let parametric: i32 = reader.parse(block_header[2], "parametric flag")?;
        let count: usize = reader.parse(block_header[3], "node count")?;

        if parametric != 0 {
            return Err(MeshError::UnsupportedParametricNodes { line: reader.line });
        }
        if count == 0 {
            continue;
        }

        // Auto-detect: read the first data line to determine the format.
        let first = reader.read_line()?;
        let first_tokens: Vec<&str> = first.split_whitespace().collect();
        if first_tokens.is_empty() {
            return Err(reader.malformed("empty node line"));
        }

        // If the first line has one token, it is the separated format: all
        // tags come first, then all coordinates. Otherwise it is the compact
        // format: "tag x y z" per line.
        if first_tokens.len() == 1 {
            let mut tags = Vec::with_capacity(count);
            tags.push(reader.parse::<usize>(first_tokens[0], "node tag")?);

            while tags.len() < count {
                let line = reader.read_line()?;
                for token in line.split_whitespace() {
                    if tags.len() == count {
                        break;
                    }
                    tags.push(reader.parse::<usize>(token, "node tag")?);
                }
            }

            for tag in tags {
                let line = reader.read_line()?;
                let tokens: Vec<&str> = line.split_whitespace().collect();
                if tokens.len() < D {
                    return Err(reader.malformed(format!("expected at least {D} node coordinates")));
                }

                let mut coords = [0.0_f64; D];
                for (d, token) in tokens.iter().take(D).enumerate() {
                    coords[d] = reader.parse(token, "node coordinate")?;
                }
                insert_node(&mut nodes, tag, coords, reader.line)?;
            }
        } else {
            parse_compact_node_line::<D>(reader, &first, &mut nodes)?;
            for _ in 1..count {
                let line = reader.read_line()?;
                parse_compact_node_line::<D>(reader, &line, &mut nodes)?;
            }
        }
    }

    reader.expect_end("$EndNodes")?;

    if nodes.len() != num_nodes {
        return Err(MeshError::NodeCountMismatch {
            expected: num_nodes,
            found: nodes.len(),
        });
    }

    Ok(nodes)
}

fn insert_node<const D: usize>(
    nodes: &mut HashMap<usize, Node<D>>,
    tag: usize,
    coords: [f64; D],
    line: usize,
) -> Result<(), MeshError> {
    let id = tag
        .checked_sub(1)
        .ok_or(MeshError::InvalidNodeTag { tag, line })?;

    if nodes.insert(id, Node::new(id, coords)).is_some() {
        return Err(MeshError::DuplicateNode { tag, line });
    }
    Ok(())
}

fn parse_compact_node_line<const D: usize>(
    reader: &MshReader,
    line: &str,
    nodes: &mut HashMap<usize, Node<D>>,
) -> Result<(), MeshError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 1 + D {
        return Err(reader.malformed(format!("expected tag + {D} node coordinates")));
    }

    let tag: usize = reader.parse(tokens[0], "node tag")?;
    let mut coords = [0.0_f64; D];
    for (d, token) in tokens.iter().skip(1).take(D).enumerate() {
        coords[d] = reader.parse(token, "node coordinate")?;
    }

    insert_node(nodes, tag, coords, reader.line)
}

// ============================================================================
// $Elements (v4)
// ============================================================================

/// Read raw element connectivity grouped by physical group, without constructing
/// concrete element types. Each block's element type and per-element node ids are
/// preserved so callers can materialize typed vectors on demand.
fn read_elements_v4(
    reader: &mut MshReader,
    entity_to_phys: &HashMap<(i32, i32), i32>,
) -> Result<HashMap<i32, ElementBlock>, MeshError> {
    let header: Vec<usize> = reader
        .read_line()?
        .split_whitespace()
        .map(|s| reader.parse(s, "element header"))
        .collect::<Result<_, _>>()?;

    if header.len() != 4 {
        return Err(reader.malformed("expected 4 $Elements header fields"));
    }
    let num_entity_blocks = header[0];
    let num_elements = header[1];

    let mut blocks: HashMap<i32, ElementBlock> = HashMap::new();
    let mut total = 0usize;

    for _ in 0..num_entity_blocks {
        let block_line = reader.read_line()?;
        let block_header: Vec<&str> = block_line.split_whitespace().collect();
        if block_header.len() < 4 {
            return Err(reader.malformed("malformed element entity block header"));
        }

        let entity_dim: i32 = reader.parse(block_header[0], "entityDim")?;
        let entity_tag: i32 = reader.parse(block_header[1], "entityTag")?;
        let elem_type: i32 = reader.parse(block_header[2], "elementType")?;
        let count: usize = reader.parse(block_header[3], "element count")?;

        let phys_group = entity_to_phys
            .get(&(entity_dim, entity_tag))
            .copied()
            .unwrap_or(0);

        let block = match blocks.entry(phys_group) {
            Entry::Occupied(entry) => {
                let existing = entry.into_mut();
                if existing.elem_type != elem_type {
                    return Err(MeshError::ElementTypeMismatch {
                        group: phys_group,
                        expected: existing.elem_type,
                        found: elem_type,
                    });
                }
                existing
            }
            Entry::Vacant(entry) => entry.insert(ElementBlock {
                elem_type,
                connectivity: Vec::with_capacity(count),
            }),
        };

        block.connectivity.reserve(count);

        for _ in 0..count {
            let line = reader.read_line()?;
            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.len() < 2 {
                return Err(reader.malformed("expected element tag and at least one node tag"));
            }

            let _element_tag: usize = reader.parse(tokens[0], "element tag")?;
            let mut ids = Vec::with_capacity(tokens.len() - 1);
            for token in &tokens[1..] {
                let tag: usize = reader.parse(token, "node tag")?;
                let id = tag.checked_sub(1).ok_or(MeshError::InvalidNodeTag {
                    tag,
                    line: reader.line,
                })?;
                ids.push(id);
            }

            block.connectivity.push(ids);
            total += 1;
        }
    }

    reader.expect_end("$EndElements")?;

    if total != num_elements {
        return Err(MeshError::ElementCountMismatch {
            expected: num_elements,
            found: total,
        });
    }

    for block in blocks.values_mut() {
        block.connectivity.shrink_to_fit();
    }

    Ok(blocks)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::Element;
    use std::io::Write;

    #[test]
    fn reads_segments_and_triangles() {
        let mesh = read_gmsh::<2>("examples/msh/patchtest_tri3.msh").unwrap();
        let segs: Vec<Seg2<2>> = mesh.elements("Γᵍ").unwrap();
        assert_eq!(segs.len(), 40, "expected 40 segment elements");
        let tris: Vec<Tri3<2>> = mesh.elements("Ω").unwrap();
        assert_eq!(tris.len(), 288, "expected 288 triangle elements");
    }

    #[test]
    fn element_type_reports_group_element_type() {
        let mesh = read_gmsh::<2>("examples/msh/patchtest_tri3.msh").unwrap();
        assert_eq!(
            mesh.element_type("Ω").unwrap(),
            Tri3::<2>::GMSH_ELEMENT_TYPE
        );
        assert_eq!(
            mesh.element_type("Γᵍ").unwrap(),
            Seg2::<2>::GMSH_ELEMENT_TYPE
        );
    }

    #[test]
    fn element_type_mismatch_returns_error() {
        let mesh = read_gmsh::<2>("examples/msh/patchtest_tri3.msh").unwrap();
        assert!(matches!(
            mesh.elements::<Quad4<2>>("Ω"),
            Err(MeshError::ElementTypeMismatch { .. })
        ));
    }

    #[test]
    fn unknown_physical_name_returns_error() {
        let mesh = read_gmsh::<2>("examples/msh/patchtest_tri3.msh").unwrap();
        assert!(matches!(
            mesh.elements::<Tri3<2>>("does-not-exist"),
            Err(MeshError::UnknownPhysicalName { .. })
        ));
    }

    #[test]
    fn quad4_from_gmsh() {
        let mut nodes = HashMap::new();
        nodes.insert(0, Node::new(0, [0.0, 0.0]));
        nodes.insert(1, Node::new(1, [1.0, 0.0]));
        nodes.insert(2, Node::new(2, [1.0, 1.0]));
        nodes.insert(3, Node::new(3, [0.0, 1.0]));
        let table = NodeTable::from_map(&nodes).unwrap();
        let q = Quad4::from_gmsh(&[0, 1, 2, 3], &table).unwrap();
        assert_eq!(q.id(), [0, 1, 2, 3]);
    }

    #[test]
    fn from_gmsh_missing_node_returns_error() {
        let nodes: HashMap<usize, Node<2>> = HashMap::new();
        let table = NodeTable::from_map(&nodes).unwrap();
        assert!(matches!(
            Quad4::from_gmsh(&[0, 1, 2, 3], &table),
            Err(MeshError::NodeIndexOutOfRange { .. })
        ));
    }

    #[test]
    fn reads_points_and_segments() {
        let mesh = read_gmsh::<1>("examples/msh/bar1d.msh").unwrap();
        // phys 0 (unnamed) → 2 point elements; phys 1 "line" → 1 segment element.
        let points: Vec<Poi1<1>> = mesh.elements_by_tag(0).unwrap();
        assert_eq!(points.len(), 2, "expected 2 point elements");
        let segs: Vec<Seg2<1>> = mesh.elements("line").unwrap();
        assert_eq!(segs.len(), 1, "expected 1 segment element");
    }

    #[test]
    fn rejects_non_v4_version() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "$MeshFormat").unwrap();
        writeln!(file, "2.2 0 8").unwrap();
        writeln!(file, "$EndMeshFormat").unwrap();

        assert!(matches!(
            read_gmsh::<2>(file.path()),
            Err(MeshError::UnsupportedVersion { .. })
        ));
    }

    #[test]
    fn rejects_binary_format() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "$MeshFormat").unwrap();
        writeln!(file, "4.1 1 8").unwrap();
        writeln!(file, "$EndMeshFormat").unwrap();

        assert!(matches!(
            read_gmsh::<2>(file.path()),
            Err(MeshError::BinaryNotSupported)
        ));
    }

    #[test]
    fn elements_missing_node_returns_error() {
        let mut blocks = HashMap::new();
        blocks.insert(
            1,
            ElementBlock {
                elem_type: 2,
                connectivity: vec![vec![0, 1, 2]],
            },
        );
        let mesh = GmshMesh::<2> {
            nodes: HashMap::new(),
            physical_names: HashMap::from([(1, "Ω".to_string())]),
            name_to_tag: HashMap::from([("Ω".to_string(), 1)]),
            blocks,
        };

        assert!(matches!(
            mesh.elements::<Tri3<2>>("Ω"),
            Err(MeshError::MissingNode { .. })
        ));
    }

    #[test]
    fn missing_entities_section_returns_error() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "$MeshFormat").unwrap();
        writeln!(file, "4.1 0 8").unwrap();
        writeln!(file, "$EndMeshFormat").unwrap();

        assert!(matches!(
            read_gmsh::<2>(file.path()),
            Err(MeshError::MissingSection {
                section: "$Entities"
            })
        ));
    }

    #[test]
    fn malformed_physical_names_returns_error() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "$MeshFormat").unwrap();
        writeln!(file, "4.1 0 8").unwrap();
        writeln!(file, "$EndMeshFormat").unwrap();
        writeln!(file, "$PhysicalNames").unwrap();
        writeln!(file, "1").unwrap();
        writeln!(file, "2 oops").unwrap();
        writeln!(file, "$EndPhysicalNames").unwrap();

        assert!(matches!(
            read_gmsh::<2>(file.path()),
            Err(MeshError::MalformedLine { .. })
        ));
    }

    #[test]
    fn physical_groups_are_sorted_and_counted() {
        let mesh = read_gmsh::<2>("examples/msh/patchtest_tri3.msh").unwrap();
        let groups: Vec<_> = mesh.physical_groups().collect();
        assert!(groups.windows(2).all(|pair| pair[0].tag <= pair[1].tag));

        let omega = mesh.physical_group_by_name("Ω").unwrap();
        assert_eq!(omega.name, Some("Ω"));
        assert_eq!(omega.element_type, Tri3::<2>::GMSH_ELEMENT_TYPE);
        assert_eq!(omega.element_count, 288);
        assert_eq!(mesh.physical_group(omega.tag).unwrap(), omega);
        assert!(mesh.has_physical_group("Ω"));
        assert!(!mesh.has_physical_group("does-not-exist"));

        let names: Vec<_> = mesh.physical_names().collect();
        assert!(names.windows(2).all(|pair| pair[0].0 <= pair[1].0));

        let table = mesh.node_table().unwrap();
        let tris: Vec<Tri3<2>> = mesh.elements_by_tag_with_table(omega.tag, &table).unwrap();
        assert_eq!(tris.len(), 288);
    }

    #[test]
    fn unknown_physical_tag_returns_error() {
        let mesh = read_gmsh::<2>("examples/msh/patchtest_tri3.msh").unwrap();
        assert!(matches!(
            mesh.physical_group(-1),
            Err(MeshError::UnknownPhysicalTag { .. })
        ));
    }

    #[test]
    fn duplicate_physical_names_return_error() {
        let names = HashMap::from([(1, "wall".to_string()), (2, "wall".to_string())]);
        assert!(matches!(
            build_name_to_tag(&names),
            Err(MeshError::DuplicatePhysicalName { .. })
        ));
    }

    #[test]
    fn node_table_rejects_non_contiguous_ids() {
        let nodes = HashMap::from([(1, Node::new(1, [0.0, 0.0]))]);
        assert!(matches!(
            NodeTable::from_map(&nodes),
            Err(MeshError::NodeIndexOutOfRange { .. })
        ));
    }

    #[test]
    fn node_table_get_rejects_out_of_range() {
        let table = NodeTable::from_nodes(vec![Node::new(0, [0.0, 0.0])]).unwrap();
        assert!(matches!(
            table.get(1),
            Err(MeshError::NodeIndexOutOfRange { .. })
        ));
    }

    #[test]
    fn gmsh_mesh_debug_is_a_summary() {
        let mesh = read_gmsh::<2>("examples/msh/patchtest_tri3.msh").unwrap();
        let text = format!("{mesh:?}");
        assert!(text.contains("GmshMesh"));
        assert!(text.contains("n_nodes"));
        assert!(text.contains("groups"));
    }
}
