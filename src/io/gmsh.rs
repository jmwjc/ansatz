use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use crate::element::Element;
use crate::element::cr::{CRSeg3, CRTet4, CRTri3, CRTri4};
use crate::element::poi1::Poi1;
use crate::element::quad4::Quad4;
use crate::element::seg2::Seg2;
use crate::element::tet4::Tet4;
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
    EdgeNumberingMissing,
    MissingEdge {
        node_a: usize,
        node_b: usize,
    },
    InvalidFacet {
        index: usize,
        facets: usize,
    },
    FaceNumberingMissing,
    MissingFace {
        node_a: usize,
        node_b: usize,
        node_c: usize,
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
            Self::EdgeNumberingMissing => write!(
                f,
                "the mesh has no edge numbering; read it with `GmshOptions {{ edges: true }}` or pass an explicit edge table"
            ),
            Self::MissingEdge { node_a, node_b } => write!(
                f,
                "the edge ({node_a}, {node_b}) has no degree of freedom in the edge table"
            ),
            Self::InvalidFacet { index, facets } => write!(
                f,
                "facet index {index} out of range for an element with {facets} facets"
            ),
            Self::FaceNumberingMissing => write!(
                f,
                "the mesh has no face numbering; read it with `GmshOptions {{ faces: true }}`"
            ),
            Self::MissingFace {
                node_a,
                node_b,
                node_c,
            } => write!(
                f,
                "the face ({node_a}, {node_b}, {node_c}) has no degree of freedom in the face table"
            ),
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
/// The table holds no gaps and no redundant `Node::id` field, so it can only be
/// built from node ids that are exactly `0..len`. Reading with
/// [`GmshOptions::renumber`] guarantees that; without renumbering, a mesh whose
/// GMSH tags leave gaps is rejected by [`NodeTable::from_map`].
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

/// Edge table used to number the edge degrees of freedom of a
/// Crouzeix–Raviart ([`CRTri3`]) discretization.
///
/// A CR dof belongs to an **edge** of the mesh, not to a vertex, and
/// [`CRTri3::from_tri3`] expects the global dof id of each edge. This table assigns
/// such ids on demand: the first time an edge (given by its two end vertex ids, in
/// either order) is seen it receives `offset + k`, so the edge dofs occupy the
/// contiguous block `offset..offset + len`.
///
/// The edge numbering is **independent of the node numbering**: a CR dof
/// discretizes a different quantity than any vertex-based field, so the two index
/// ranges may overlap without conflict (`GmshMesh::edge_table` and
/// [`GmshOptions::edges`] start at `0`). Pass a non-zero `offset` only when the CR
/// dofs must be placed after another field's dofs in one global system.
///
/// ```
/// use ansatz::io::gmsh::EdgeTable;
///
/// let mut edges = EdgeTable::default();
/// assert_eq!(edges.id(0, 1), 0);
/// assert_eq!(edges.id(1, 0), 0); // same edge, either orientation
/// assert_eq!(edges.id(1, 2), 1);
/// assert_eq!(edges.len(), 2);
///
/// // Block placement for a mixed (multi-field) system:
/// let mut shifted = EdgeTable::new(10);
/// assert_eq!(shifted.id(0, 1), 10);
/// ```
#[derive(Clone, Debug, Default)]
pub struct EdgeTable {
    edges: HashMap<(usize, usize), usize>,
    offset: usize,
}

impl EdgeTable {
    /// Creates an empty table whose edge dofs start at `offset`.
    ///
    /// `offset` is `0` for a pure CR discretization (`EdgeTable::default()`); use a
    /// non-zero value to place the CR dofs after another field's dofs.
    pub fn new(offset: usize) -> Self {
        EdgeTable {
            edges: HashMap::new(),
            offset,
        }
    }

    /// First global id handed out to an edge.
    pub fn offset(&self) -> usize {
        self.offset
    }

    /// Number of distinct edges seen so far.
    pub fn len(&self) -> usize {
        self.edges.len()
    }

    /// Whether no edge has been registered yet.
    pub fn is_empty(&self) -> bool {
        self.edges.is_empty()
    }

    /// Global dof id of the edge between vertices `a` and `b`, registering the
    /// edge (and thereby allocating the next id) if it is new.
    pub fn id(&mut self, a: usize, b: usize) -> usize {
        let key = if a <= b { (a, b) } else { (b, a) };
        let next = self.offset + self.edges.len();
        *self.edges.entry(key).or_insert(next)
    }

    /// Global dof id of an already registered edge, or `None`.
    pub fn get(&self, a: usize, b: usize) -> Option<usize> {
        let key = if a <= b { (a, b) } else { (b, a) };
        self.edges.get(&key).copied()
    }
}

/// Face table used to number the face degrees of freedom of a
/// Crouzeix–Raviart [`CRTet4`] discretization.
///
/// The 3D counterpart of [`EdgeTable`]: a CR tetrahedron has one dof per face, and
/// [`CRTet4::from_tet4`] expects the global id of each face. Ids are assigned on
/// demand to a face given by its three vertex ids (in any order), so the dofs occupy
/// `offset..offset + len`.
///
/// Like the edge numbering, the face numbering is independent of the node numbering.
///
/// ```
/// use ansatz::io::gmsh::FaceTable;
///
/// let mut faces = FaceTable::default();
/// assert_eq!(faces.id(0, 1, 2), 0);
/// assert_eq!(faces.id(2, 0, 1), 0); // same face, any order
/// assert_eq!(faces.id(1, 2, 3), 1);
/// assert_eq!(faces.len(), 2);
/// ```
#[derive(Clone, Debug, Default)]
pub struct FaceTable {
    faces: HashMap<[usize; 3], usize>,
    offset: usize,
}

impl FaceTable {
    /// Creates an empty table whose face dofs start at `offset`.
    ///
    /// `offset` is `0` for a pure CR discretization; use a non-zero value to place
    /// the CR dofs after another field's dofs.
    pub fn new(offset: usize) -> Self {
        FaceTable {
            faces: HashMap::new(),
            offset,
        }
    }

    /// First global id handed out to a face.
    pub fn offset(&self) -> usize {
        self.offset
    }

    /// Number of distinct faces seen so far.
    pub fn len(&self) -> usize {
        self.faces.len()
    }

    /// Whether no face has been registered yet.
    pub fn is_empty(&self) -> bool {
        self.faces.is_empty()
    }

    /// Global dof id of the face through vertices `a`, `b` and `c`, registering the
    /// face (and thereby allocating the next id) if it is new.
    pub fn id(&mut self, a: usize, b: usize, c: usize) -> usize {
        let key = canonical_face([a, b, c]);
        let next = self.offset + self.faces.len();
        *self.faces.entry(key).or_insert(next)
    }

    /// Global dof id of an already registered face, or `None`.
    pub fn get(&self, a: usize, b: usize, c: usize) -> Option<usize> {
        self.faces.get(&canonical_face([a, b, c])).copied()
    }
}

/// Sorts a face's three vertex ids so that a face has a single key.
fn canonical_face(mut ids: [usize; 3]) -> [usize; 3] {
    ids.sort_unstable();
    ids
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
    /// Edge numbering created by [`GmshOptions::edges`], if requested.
    edges: Option<EdgeTable>,
    /// Parent element of every triangle edge, created by [`GmshOptions::edges`].
    parents: Option<HashMap<(usize, usize), FacetParent>>,
    /// Face numbering created by [`GmshOptions::faces`], if requested.
    faces: Option<FaceTable>,
    /// Parent tetrahedron of every tetrahedron face, created by [`GmshOptions::faces`].
    face_parents: Option<HashMap<[usize; 3], FaceParent>>,
}

/// Parent tetrahedron of a mesh face: the CR trace element is built from it.
#[derive(Clone, Copy, Debug)]
struct FaceParent {
    /// The parent's four face dofs, in `CRTet4::id()` order.
    facets: [usize; 4],
    /// The parent's four vertex ids, in local order.
    vertices: [usize; 4],
    /// Local facet index: the facet opposite parent vertex `facet`.
    facet: usize,
}

/// Parent triangle of a mesh edge: the CR trace element is built from it.
#[derive(Clone, Copy, Debug)]
struct FacetParent {
    /// The parent's three facet dofs, in `CRTri3::id()` order.
    facets: [usize; 3],
    /// The parent's three vertex ids, in local order.
    vertices: [usize; 3],
    /// Local facet index: the facet opposite parent vertex `facet`.
    facet: usize,
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
            .field("n_edges", &self.edges.as_ref().map(EdgeTable::len))
            .field("n_faces", &self.faces.as_ref().map(FaceTable::len))
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

impl<const D: usize> GmshMesh<D> {
    /// Creates an empty [`EdgeTable`] numbering edges from `0`.
    ///
    /// The edge numbering is independent of the node numbering (see [`EdgeTable`]),
    /// so no offset is needed; use [`EdgeTable::new`] instead if the CR dofs must
    /// follow another field's dofs in a mixed global system. Pass the table to the
    /// `crtri3_elements_with_table*` readers. If the mesh was read with
    /// [`GmshOptions::edges`], the numbering already exists — see
    /// [`GmshMesh::edge_numbering`].
    pub fn edge_table(&self) -> EdgeTable {
        EdgeTable::default()
    }

    /// The edge numbering created by [`GmshOptions::edges`], if it was requested.
    ///
    /// `offset() .. offset() + len()` is the range of CR edge dofs; for a pure CR
    /// discretization that is `0 .. n_edges`, i.e. the CR unknown count.
    pub fn edge_numbering(&self) -> Option<&EdgeTable> {
        self.edges.as_ref()
    }

    /// Builds the Crouzeix–Raviart elements of a physical name.
    ///
    /// This is the CR counterpart of [`GmshMesh::elements`] and requires the mesh
    /// to carry an edge numbering, i.e. to be read with
    /// [`GmshOptions::edges`]. The group must hold GMSH triangles (type 2).
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::UnknownPhysicalName`] if the name is not present,
    /// [`MeshError::EmptyPhysicalGroup`] if that group has no elements,
    /// [`MeshError::ElementTypeMismatch`] if the group holds another element type,
    /// [`MeshError::EdgeNumberingMissing`] if the mesh was read without
    /// [`GmshOptions::edges`], or a connectivity error if an element cannot be
    /// built.
    pub fn crtri3_elements(&self, physical_name: &str) -> Result<Vec<CRTri3<D>>, MeshError> {
        let tag = self.physical_tag(physical_name)?;
        self.crtri3_elements_by_tag(tag)
    }

    /// Builds the Crouzeix–Raviart elements of a physical tag.
    ///
    /// # Errors
    ///
    /// Same as [`GmshMesh::crtri3_elements`].
    pub fn crtri3_elements_by_tag(&self, tag: i32) -> Result<Vec<CRTri3<D>>, MeshError> {
        let edges = self.edges.as_ref().ok_or(MeshError::EdgeNumberingMissing)?;
        let node_table = self.node_table()?;
        let triangles: Vec<Tri3<D>> = self.elements_by_tag_with_table(tag, &node_table)?;

        triangles
            .iter()
            .map(|tri| {
                CRTri3::try_from_tri3(tri, |a, b| {
                    edges.get(a, b).ok_or(MeshError::MissingEdge {
                        node_a: a,
                        node_b: b,
                    })
                })
            })
            .collect()
    }

    /// Builds the Crouzeix–Raviart **trace** elements of a boundary group.
    ///
    /// Boundary conditions for a CR discretization use [`CRSeg3`], the restriction
    /// of the CR basis to one facet of a triangle: its dofs are the **three** edge
    /// dofs of the adjacent triangle, because the field on a facet depends on all of
    /// them (a plain [`Seg2`] would carry vertex values, which CR does not have).
    /// The reader therefore groups the boundary segments (GMSH type 1) with the
    /// triangle they belong to and derives the local facet index, so the resulting
    /// element can be passed straight to [`crate::operation::heat::heat_flux`],
    /// [`crate::operation::heat::convection`],
    /// [`crate::operation::heat::temperature_penalty`],
    /// [`crate::operation::heat::temperature_lagrangian`],
    /// [`crate::operation::elasticity::traction`], … like any other boundary
    /// element.
    ///
    /// Requires a mesh read with [`GmshOptions::edges`] (which records the parent
    /// triangle of every edge).
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::UnknownPhysicalName`] if the name is not present,
    /// [`MeshError::EmptyPhysicalGroup`] if that group has no elements,
    /// [`MeshError::ElementTypeMismatch`] if the group does not hold segments,
    /// [`MeshError::EdgeNumberingMissing`] if the mesh was read without
    /// [`GmshOptions::edges`], [`MeshError::MissingEdge`] if a segment is not an edge
    /// of any triangle of the mesh, or a connectivity error if an element cannot be
    /// built.
    pub fn crtri3_boundary(&self, physical_name: &str) -> Result<Vec<CRSeg3<D>>, MeshError> {
        let tag = self.physical_tag(physical_name)?;
        self.crtri3_boundary_by_tag(tag)
    }

    /// Builds the Crouzeix–Raviart trace elements of a physical tag.
    ///
    /// # Errors
    ///
    /// Same as [`GmshMesh::crtri3_boundary`].
    pub fn crtri3_boundary_by_tag(&self, tag: i32) -> Result<Vec<CRSeg3<D>>, MeshError> {
        let parents = self
            .parents
            .as_ref()
            .ok_or(MeshError::EdgeNumberingMissing)?;
        let node_table = self.node_table()?;
        let segments: Vec<Seg2<D>> = self.elements_by_tag_with_table(tag, &node_table)?;

        let mut facets = Vec::with_capacity(segments.len());
        for segment in &segments {
            let [a, b] = segment.id();
            let key = (a.min(b), a.max(b));
            let parent = parents.get(&key).ok_or(MeshError::MissingEdge {
                node_a: key.0,
                node_b: key.1,
            })?;

            let mut vertices = [Node::new(0, [0.0; D]); 3];
            for (index, node_id) in parent.vertices.iter().enumerate() {
                vertices[index] = Node::new(*node_id, node_table.get(*node_id)?);
            }

            facets.push(CRSeg3::new(parent.facets, vertices, parent.facet).ok_or(
                MeshError::InvalidFacet {
                    index: parent.facet,
                    facets: 3,
                },
            )?);
        }

        Ok(facets)
    }

    /// Creates an empty [`FaceTable`] numbering faces from `0`.
    ///
    /// Use it when you want to number faces yourself; the 3D CR reader
    /// [`GmshMesh::crtet4_elements`] draws its numbering from
    /// [`GmshOptions::faces`] instead.
    pub fn face_table(&self) -> FaceTable {
        FaceTable::default()
    }

    /// The face numbering created by [`GmshOptions::faces`], if it was requested.
    ///
    /// `offset() .. offset() + len()` is the range of CR face dofs; for a pure CR
    /// discretization that is `0 .. n_faces`.
    pub fn face_numbering(&self) -> Option<&FaceTable> {
        self.faces.as_ref()
    }

    /// Builds the Crouzeix–Raviart tetrahedra of a physical name.
    ///
    /// The 3D counterpart of [`GmshMesh::crtri3_elements`]: the group must hold GMSH
    /// tetrahedra (type 4), which are read as [`Tet4`] and then converted with
    /// [`CRTet4::from_tet4`]; the face dofs come from the numbering created by
    /// [`GmshOptions::faces`].
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::UnknownPhysicalName`] if the name is not present,
    /// [`MeshError::EmptyPhysicalGroup`] if that group has no elements,
    /// [`MeshError::ElementTypeMismatch`] if the group holds another element type,
    /// [`MeshError::FaceNumberingMissing`] if the mesh was read without
    /// [`GmshOptions::faces`], or a connectivity error.
    pub fn crtet4_elements(&self, physical_name: &str) -> Result<Vec<CRTet4<D>>, MeshError> {
        let tag = self.physical_tag(physical_name)?;
        self.crtet4_elements_by_tag(tag)
    }

    /// Builds the Crouzeix–Raviart tetrahedra of a physical tag.
    ///
    /// # Errors
    ///
    /// Same as [`GmshMesh::crtet4_elements`].
    pub fn crtet4_elements_by_tag(&self, tag: i32) -> Result<Vec<CRTet4<D>>, MeshError> {
        let faces = self.faces.as_ref().ok_or(MeshError::FaceNumberingMissing)?;
        let node_table = self.node_table()?;
        let tetrahedra: Vec<Tet4<D>> = self.elements_by_tag_with_table(tag, &node_table)?;

        tetrahedra
            .iter()
            .map(|tet| {
                CRTet4::try_from_tet4(tet, |face| {
                    faces
                        .get(face[0], face[1], face[2])
                        .ok_or(MeshError::MissingFace {
                            node_a: face[0],
                            node_b: face[1],
                            node_c: face[2],
                        })
                })
            })
            .collect()
    }

    /// Builds the Crouzeix–Raviart **trace** elements of a 3D boundary group.
    ///
    /// The 3D counterpart of [`GmshMesh::crtri3_boundary`]: boundary conditions for a
    /// CR tetrahedral discretization use [`CRTri4`], the restriction of the CR basis
    /// to one face of a tetrahedron, whose dofs are the **four** face dofs of the
    /// adjacent tetrahedron. The reader groups the boundary triangles (GMSH type 2)
    /// with the tetrahedron they belong to and derives the local facet index.
    ///
    /// Requires a mesh read with [`GmshOptions::faces`].
    ///
    /// # Errors
    ///
    /// Returns [`MeshError::UnknownPhysicalName`] if the name is not present,
    /// [`MeshError::EmptyPhysicalGroup`] if that group has no elements,
    /// [`MeshError::ElementTypeMismatch`] if the group does not hold triangles,
    /// [`MeshError::FaceNumberingMissing`] if the mesh was read without
    /// [`GmshOptions::faces`], or [`MeshError::MissingFace`] if a triangle is not a
    /// face of any tetrahedron of the mesh.
    pub fn crtet4_boundary(&self, physical_name: &str) -> Result<Vec<CRTri4<D>>, MeshError> {
        let tag = self.physical_tag(physical_name)?;
        self.crtet4_boundary_by_tag(tag)
    }

    /// Builds the Crouzeix–Raviart trace elements of a physical tag.
    ///
    /// # Errors
    ///
    /// Same as [`GmshMesh::crtet4_boundary`].
    pub fn crtet4_boundary_by_tag(&self, tag: i32) -> Result<Vec<CRTri4<D>>, MeshError> {
        let face_parents = self
            .face_parents
            .as_ref()
            .ok_or(MeshError::FaceNumberingMissing)?;
        let node_table = self.node_table()?;
        let triangles: Vec<Tri3<D>> = self.elements_by_tag_with_table(tag, &node_table)?;

        let mut facets = Vec::with_capacity(triangles.len());
        for triangle in &triangles {
            let key = canonical_face(triangle.id());
            let parent = face_parents.get(&key).ok_or(MeshError::MissingFace {
                node_a: key[0],
                node_b: key[1],
                node_c: key[2],
            })?;

            let mut vertices = [Node::new(0, [0.0; D]); 4];
            for (index, node_id) in parent.vertices.iter().enumerate() {
                vertices[index] = Node::new(*node_id, node_table.get(*node_id)?);
            }

            facets.push(CRTri4::new(parent.facets, vertices, parent.facet).ok_or(
                MeshError::InvalidFacet {
                    index: parent.facet,
                    facets: 4,
                },
            )?);
        }

        Ok(facets)
    }

    /// Builds the Crouzeix–Raviart elements of a physical name using a
    /// caller-provided dense [`NodeTable`] and [`EdgeTable`].
    ///
    /// Unlike [`GmshMesh::crtri3_elements`], the edge numbering comes from `edges`
    /// ([`GmshMesh::edge_table`] creates an empty one): new edges are appended in
    /// visiting order. Share one table across all groups so that an edge on a
    /// shared interface gets a single dof.
    ///
    /// # Errors
    ///
    /// Same as [`GmshMesh::crtri3_elements`].
    pub fn crtri3_elements_with_table(
        &self,
        physical_name: &str,
        node_table: &NodeTable<D>,
        edges: &mut EdgeTable,
    ) -> Result<Vec<CRTri3<D>>, MeshError> {
        let tag = self.physical_tag(physical_name)?;
        self.crtri3_elements_by_tag_with_table(tag, node_table, edges)
    }

    /// Builds the Crouzeix–Raviart elements of a physical tag using a
    /// caller-provided dense [`NodeTable`] and [`EdgeTable`].
    ///
    /// # Errors
    ///
    /// Same as [`GmshMesh::crtri3_elements_with_table`].
    pub fn crtri3_elements_by_tag_with_table(
        &self,
        tag: i32,
        node_table: &NodeTable<D>,
        edges: &mut EdgeTable,
    ) -> Result<Vec<CRTri3<D>>, MeshError> {
        let triangles: Vec<Tri3<D>> = self.elements_by_tag_with_table(tag, node_table)?;
        Ok(triangles
            .iter()
            .map(|tri| CRTri3::from_tri3(tri, |a, b| edges.id(a, b)))
            .collect())
    }
}

/// Numbers the faces of every tetrahedron block from `0` and records the parent
/// tetrahedron of each face (needed to build CR trace elements).
///
/// The 3D counterpart of [`number_edges`]: blocks are visited in physical-tag order
/// and elements in file order, so the numbering is deterministic.
fn number_faces<const D: usize>(
    blocks: &HashMap<i32, ElementBlock>,
) -> (FaceTable, HashMap<[usize; 3], FaceParent>) {
    let mut faces = FaceTable::default();
    let mut parents: HashMap<[usize; 3], FaceParent> = HashMap::new();

    let mut tags: Vec<i32> = blocks.keys().copied().collect();
    tags.sort_unstable();
    for tag in tags {
        let Some(block) = blocks.get(&tag) else {
            continue;
        };
        if block.elem_type != Tet4::<D>::GMSH_ELEMENT_TYPE {
            continue;
        }
        for ids in &block.connectivity {
            if ids.len() < 4 {
                continue;
            }
            let vertices = [ids[0], ids[1], ids[2], ids[3]];
            // Facets in `CRTet4::id()` order: opposite vertex 0, 1, 2, 3.
            let local_faces = [(1_usize, 2_usize, 3_usize), (0, 2, 3), (0, 1, 3), (0, 1, 2)];
            let facets = std::array::from_fn(|facet| {
                let (a, b, c) = local_faces[facet];
                faces.id(vertices[a], vertices[b], vertices[c])
            });
            for (facet, (a, b, c)) in local_faces.into_iter().enumerate() {
                let key = canonical_face([vertices[a], vertices[b], vertices[c]]);
                parents.entry(key).or_insert(FaceParent {
                    facets,
                    vertices,
                    facet,
                });
            }
        }
    }

    (faces, parents)
}

/// Options controlling the preprocessing steps applied by [`read_gmsh_with`].
///
/// The default ([`GmshOptions::default`]) applies no preprocessing: the node ids
/// are the raw GMSH ones (`tag - 1`) and no edge numbering is built, exactly like
/// [`read_gmsh`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GmshOptions {
    /// Renumber node ids densely to `0..n_nodes`, in ascending raw-id order, and
    /// remap all element connectivity accordingly.
    ///
    /// GMSH node tags need not be contiguous (tags may have gaps after editing a
    /// mesh), in which case the raw ids are gapped too and
    /// [`NodeTable::from_map`] rejects them. Renumbering makes the numbering
    /// gap-free, which is what downstream preprocessing relies on — for example
    /// allocating Lagrange-multiplier dofs *after* the primary field dofs
    /// (`n_nodes + k`; see [`crate::operation::heat::temperature_lagrangian`] and
    /// [`crate::operation::elasticity::displacement_lagrangian`]).
    ///
    /// With renumbering enabled, connectivity that references a node tag absent
    /// from `$Nodes` is reported eagerly as [`MeshError::MissingNode`]; without
    /// it, the same problem surfaces later, when elements are materialized.
    pub renumber: bool,

    /// Number the edges of every triangle while reading, so that
    /// Crouzeix–Raviart elements can be materialized directly with
    /// [`GmshMesh::crtri3_elements`] (no explicit [`EdgeTable`] argument).
    ///
    /// A CR degree of freedom belongs to an edge, so the mesh needs an edge
    /// numbering in addition to the node numbering. Enabling this option creates
    /// one at read time, spanning all triangle groups (so an interior edge shared
    /// by two subdomains gets a single dof), and exposes it through
    /// [`GmshMesh::edge_numbering`]. The edge ids are `0 .. n_edges`; they are
    /// independent of the node ids, since edges and vertices index different
    /// unknowns.
    ///
    /// Note that materializing *any* typed element (including CR) still needs a
    /// dense node id range, which is why `renumber: true` is normally passed
    /// alongside; a gapped mesh reads fine but cannot be materialized.
    pub edges: bool,

    /// Number the faces of every tetrahedron while reading, so that 3D
    /// Crouzeix–Raviart elements can be materialized with
    /// [`GmshMesh::crtet4_elements`] and their boundary traces with
    /// [`GmshMesh::crtet4_boundary`].
    ///
    /// This is the 3D counterpart of [`GmshOptions::edges`]: face ids are independent
    /// of the node and edge numbering, span all tetrahedron groups (so an interior
    /// face shared by two tetrahedra gets a single dof) and are exposed through
    /// [`GmshMesh::face_numbering`].
    pub faces: bool,
}

/// Parses a GMSH v4 `.msh` file and returns the intermediate product [`GmshMesh`].
///
/// This function only reads the file and does not construct concrete element
/// types; element vectors are materialized on demand by [`GmshMesh::elements`]
/// per physical name.
///
/// Node ids are the raw GMSH ones (`tag - 1`) and no edge numbering is built; use
/// [`read_gmsh_with`] to enable the preprocessing steps of [`GmshOptions`] (dense
/// node renumbering, edge numbering for [Crouzeix–Raviart](crate::element::cr::CRTri3)
/// elements).
///
/// # Errors
///
/// Returns a [`MeshError`] if the file cannot be opened, is not GMSH v4 ASCII,
/// is malformed, or lacks a required section.
pub fn read_gmsh<const D: usize>(filepath: impl AsRef<Path>) -> Result<GmshMesh<D>, MeshError> {
    read_gmsh_with(filepath, GmshOptions::default())
}

/// Parses a GMSH v4 `.msh` file with explicit preprocessing [`GmshOptions`].
///
/// This is [`read_gmsh`] with the option to renumber the nodes densely (see
/// [`GmshOptions::renumber`]); callers that need a gap-free node numbering, for
/// instance to append Lagrange-multiplier dofs, should pass
/// `GmshOptions { renumber: true, edges: false, faces: false }`.
///
/// # Errors
///
/// Returns a [`MeshError`] if the file cannot be opened, is not GMSH v4 ASCII,
/// is malformed, or lacks a required section. With
/// [`GmshOptions::renumber`] enabled it additionally returns
/// [`MeshError::MissingNode`] if an element references a node tag that is absent
/// from `$Nodes`.
pub fn read_gmsh_with<const D: usize>(
    filepath: impl AsRef<Path>,
    options: GmshOptions,
) -> Result<GmshMesh<D>, MeshError> {
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

    // ---------- preprocessing ----------
    let (nodes, blocks) = if options.renumber {
        renumber_nodes(nodes, blocks)?
    } else {
        (nodes, blocks)
    };

    let (edges, parents) = if options.edges {
        let (edges, parents) = number_edges::<D>(&blocks);
        (Some(edges), Some(parents))
    } else {
        (None, None)
    };

    let (faces, face_parents) = if options.faces {
        let (faces, face_parents) = number_faces::<D>(&blocks);
        (Some(faces), Some(face_parents))
    } else {
        (None, None)
    };

    Ok(GmshMesh {
        nodes,
        physical_names,
        name_to_tag,
        blocks,
        edges,
        parents,
        faces,
        face_parents,
    })
}

/// Numbers the edges of every triangle block from `0` and records the parent
/// triangle of each edge (needed to build CR trace elements).
///
/// The edge numbering is independent of the node numbering (they index different
/// unknowns), so no offset is applied. Blocks are visited in physical-tag order
/// and elements in file order, so the numbering is deterministic. Groups that are
/// not triangles contribute no edges (CR boundary condition sets are usually
/// points or segments).
fn number_edges<const D: usize>(
    blocks: &HashMap<i32, ElementBlock>,
) -> (EdgeTable, HashMap<(usize, usize), FacetParent>) {
    let mut edges = EdgeTable::default();
    let mut parents: HashMap<(usize, usize), FacetParent> = HashMap::new();

    let mut tags: Vec<i32> = blocks.keys().copied().collect();
    tags.sort_unstable();
    for tag in tags {
        let Some(block) = blocks.get(&tag) else {
            continue;
        };
        if block.elem_type != Tri3::<D>::GMSH_ELEMENT_TYPE {
            continue;
        }
        for ids in &block.connectivity {
            if ids.len() < 3 {
                continue;
            }
            let vertices = [ids[0], ids[1], ids[2]];
            // Facets in `CRTri3::id()` order: opposite vertex 0, 1, 2.
            let local_facets = [(1_usize, 2_usize), (0, 2), (0, 1)];
            let facets = std::array::from_fn(|facet| {
                let (a, b) = local_facets[facet];
                edges.id(vertices[a], vertices[b])
            });
            for (facet, (a, b)) in local_facets.into_iter().enumerate() {
                let key = (vertices[a].min(vertices[b]), vertices[a].max(vertices[b]));
                parents.entry(key).or_insert(FacetParent {
                    facets,
                    vertices,
                    facet,
                });
            }
        }
    }

    (edges, parents)
}

/// `(nodes, element blocks)` pair as stored inside [`GmshMesh`].
type MeshParts<const D: usize> = (HashMap<usize, Node<D>>, HashMap<i32, ElementBlock>);

/// Renumbers node ids densely to `0..n_nodes` (ascending raw-id order) and
/// remaps every element connectivity entry.
///
/// # Errors
///
/// Returns [`MeshError::MissingNode`] if connectivity references a node id that
/// is not in `$Nodes`.
fn renumber_nodes<const D: usize>(
    nodes: HashMap<usize, Node<D>>,
    mut blocks: HashMap<i32, ElementBlock>,
) -> Result<MeshParts<D>, MeshError> {
    let mut raw_ids: Vec<usize> = nodes.keys().copied().collect();
    raw_ids.sort_unstable();

    let mut new_ids: HashMap<usize, usize> = HashMap::with_capacity(raw_ids.len());
    let mut renumbered: HashMap<usize, Node<D>> = HashMap::with_capacity(raw_ids.len());
    for (new_id, raw_id) in raw_ids.into_iter().enumerate() {
        if let Some(node) = nodes.get(&raw_id) {
            renumbered.insert(new_id, Node::new(new_id, node.coordinates));
        }
        new_ids.insert(raw_id, new_id);
    }

    // Visit the physical groups in tag order so that a shared problem reports
    // the same error every time.
    let mut tags: Vec<i32> = blocks.keys().copied().collect();
    tags.sort_unstable();
    for tag in tags {
        if let Some(block) = blocks.get_mut(&tag) {
            for (element, ids) in block.connectivity.iter_mut().enumerate() {
                for id in ids.iter_mut() {
                    *id = *new_ids.get(id).ok_or(MeshError::MissingNode {
                        group: tag,
                        element,
                        tag: *id,
                    })?;
                }
            }
        }
    }

    Ok((renumbered, blocks))
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
            edges: None,
            parents: None,
            faces: None,
            face_parents: None,
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

    /// Writes a 1D mesh whose node tags are `10, 20, 30` (so the raw ids
    /// `tag - 1` are the gapped `9, 19, 29`).
    ///
    /// It has two point groups (`$PhysicalNames` tag 2 "end", i.e. the two
    /// `Poi1` elements) and one curve group (tag 1 "line", two `Seg2` elements).
    fn write_gapped_tag_mesh() -> tempfile::NamedTempFile {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "$MeshFormat").unwrap();
        writeln!(file, "4.1 0 8").unwrap();
        writeln!(file, "$EndMeshFormat").unwrap();
        writeln!(file, "$PhysicalNames").unwrap();
        writeln!(file, "2").unwrap();
        writeln!(file, "0 2 \"end\"").unwrap();
        writeln!(file, "1 1 \"line\"").unwrap();
        writeln!(file, "$EndPhysicalNames").unwrap();
        writeln!(file, "$Entities").unwrap();
        writeln!(file, "2 1 0 0").unwrap();
        writeln!(file, "10 0 0 0 1 2").unwrap();
        writeln!(file, "30 1 0 0 1 2").unwrap();
        writeln!(file, "1 0 0 0 1 0 0 1 1 10 30").unwrap();
        writeln!(file, "$EndEntities").unwrap();
        writeln!(file, "$Nodes").unwrap();
        writeln!(file, "1 3 10 30").unwrap();
        writeln!(file, "1 1 0 3").unwrap();
        writeln!(file, "10").unwrap();
        writeln!(file, "20").unwrap();
        writeln!(file, "30").unwrap();
        writeln!(file, "0 0 0").unwrap();
        writeln!(file, "0.5 0 0").unwrap();
        writeln!(file, "1 0 0").unwrap();
        writeln!(file, "$EndNodes").unwrap();
        writeln!(file, "$Elements").unwrap();
        writeln!(file, "3 4 1 4").unwrap();
        writeln!(file, "0 10 15 1").unwrap();
        writeln!(file, "1 10").unwrap();
        writeln!(file, "0 30 15 1").unwrap();
        writeln!(file, "2 30").unwrap();
        writeln!(file, "1 1 1 2").unwrap();
        writeln!(file, "3 10 20").unwrap();
        writeln!(file, "4 20 30").unwrap();
        writeln!(file, "$EndElements").unwrap();
        file
    }

    #[test]
    fn renumbering_is_off_by_default() {
        // The default keeps the raw GMSH ids, which are gapped here, so the dense
        // node table cannot be built.
        let file = write_gapped_tag_mesh();
        let mesh = read_gmsh::<1>(file.path()).unwrap();
        let mut ids: Vec<usize> = mesh.nodes.keys().copied().collect();
        ids.sort_unstable();
        assert_eq!(ids, vec![9, 19, 29]);
        assert!(matches!(
            mesh.node_table(),
            Err(MeshError::NodeIndexOutOfRange { .. })
        ));
        assert!(matches!(
            mesh.elements::<Seg2<1>>("line"),
            Err(MeshError::MissingNode { .. } | MeshError::NodeIndexOutOfRange { .. })
        ));
        assert_eq!(
            GmshOptions::default(),
            GmshOptions {
                renumber: false,
                edges: false,
                faces: false
            }
        );
    }

    #[test]
    fn renumbering_makes_ids_dense_and_remaps_connectivity() {
        let file = write_gapped_tag_mesh();
        let mesh = read_gmsh_with::<1>(
            file.path(),
            GmshOptions {
                renumber: true,
                edges: false,
                faces: false,
            },
        )
        .unwrap();

        let table = mesh.node_table().unwrap();
        assert_eq!(table.len(), 3);
        assert_eq!(table.coords(), &[[0.0], [0.5], [1.0]]);

        // Connectivity is remapped in ascending raw-id order: ids 9, 19, 29
        // become 0, 1, 2, so the boundary points are 0 and 2 and the two
        // segments are (0, 1) and (1, 2).
        let ends: Vec<Poi1<1>> = mesh.elements("end").unwrap();
        assert_eq!(ends.len(), 2);
        assert_eq!(ends[0].id(), [0]);
        assert_eq!(ends[1].id(), [2]);

        let segs: Vec<Seg2<1>> = mesh.elements("line").unwrap();
        assert_eq!(segs.len(), 2);
        assert_eq!(segs[0].id(), [0, 1]);
        assert_eq!(segs[1].id(), [1, 2]);

        // All ids now lie in `0..n_nodes`, so multiplier dofs can be appended at
        // `n_nodes..` without colliding with mesh nodes.
        let n = table.len();
        for elm in &segs {
            assert!(elm.id().iter().all(|&id| id < n));
        }
    }

    #[test]
    fn renumbering_leaves_contiguous_meshes_unchanged() {
        let path = "examples/msh/patchtest_tri3.msh";
        let plain = read_gmsh::<2>(path).unwrap().node_table().unwrap();
        let renumbered = read_gmsh_with::<2>(
            path,
            GmshOptions {
                renumber: true,
                edges: false,
                faces: false,
            },
        )
        .unwrap()
        .node_table()
        .unwrap();
        assert_eq!(plain.coords(), renumbered.coords());
    }

    #[test]
    fn renumbering_rejects_connectivity_to_missing_nodes() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "$MeshFormat").unwrap();
        writeln!(file, "4.1 0 8").unwrap();
        writeln!(file, "$EndMeshFormat").unwrap();
        writeln!(file, "$Entities").unwrap();
        writeln!(file, "0 1 0 0").unwrap();
        writeln!(file, "1 0 0 0 1 0 0 0 0").unwrap();
        writeln!(file, "$EndEntities").unwrap();
        writeln!(file, "$Nodes").unwrap();
        writeln!(file, "1 2 1 2").unwrap();
        writeln!(file, "1 1 0 2").unwrap();
        writeln!(file, "1").unwrap();
        writeln!(file, "2").unwrap();
        writeln!(file, "0 0 0").unwrap();
        writeln!(file, "1 0 0").unwrap();
        writeln!(file, "$EndNodes").unwrap();
        writeln!(file, "$Elements").unwrap();
        writeln!(file, "1 1 1 1").unwrap();
        writeln!(file, "1 1 1 1").unwrap();
        writeln!(file, "1 1 99").unwrap();
        writeln!(file, "$EndElements").unwrap();

        assert!(matches!(
            read_gmsh_with::<1>(
                file.path(),
                GmshOptions {
                    renumber: true,
                    edges: false,
                    faces: false
                }
            ),
            Err(MeshError::MissingNode {
                group: 0,
                element: 0,
                tag: 98
            })
        ));
    }

    /// Writes a 2D mesh of two triangles sharing the diagonal:
    /// `(0,0), (1,0), (1,1), (0,1)` split by the edge between vertices 1 and 3.
    ///
    /// `tags` are the four node tags (identity is the usual `[1, 2, 3, 4]`; pass
    /// e.g. `[1, 3, 5, 7]` to obtain a mesh whose raw ids (`tag - 1`) have gaps).
    ///
    /// Groups: `Ω` (the two triangles), `Γ` (the four outer edges, in order
    /// `(1,2), (2,3), (3,4), (4,1)`), and `Γbad` (one segment from node 1 to a fifth,
    /// isolated node — not an edge of any triangle).
    fn write_two_triangle_mesh_with_tags(tags: [usize; 4]) -> tempfile::NamedTempFile {
        let [n1, n2, n3, n4] = tags;
        let stray = 90; // node tag of the isolated node
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "$MeshFormat").unwrap();
        writeln!(file, "4.1 0 8").unwrap();
        writeln!(file, "$EndMeshFormat").unwrap();
        writeln!(file, "$PhysicalNames").unwrap();
        writeln!(file, "3").unwrap();
        writeln!(file, "2 1 \"Ω\"").unwrap();
        writeln!(file, "1 2 \"Γ\"").unwrap();
        writeln!(file, "1 3 \"Γbad\"").unwrap();
        writeln!(file, "$EndPhysicalNames").unwrap();
        writeln!(file, "$Entities").unwrap();
        writeln!(file, "0 5 1 0").unwrap();
        // Five curves: the four outer edges (physical tag 2) and the stray one (3).
        writeln!(file, "1 0 0 0 1 0 0 1 2 0").unwrap();
        writeln!(file, "2 1 0 0 1 1 0 1 2 0").unwrap();
        writeln!(file, "3 0 1 0 1 1 0 1 2 0").unwrap();
        writeln!(file, "4 0 0 0 0 1 0 1 2 0").unwrap();
        writeln!(file, "5 0 0 0 0.5 2 0 1 3 0").unwrap();
        // One surface (physical tag 1).
        writeln!(file, "1 0 0 0 1 1 0 1 1 0").unwrap();
        writeln!(file, "$EndEntities").unwrap();
        writeln!(file, "$Nodes").unwrap();
        writeln!(file, "1 5 {n1} {stray}").unwrap();
        writeln!(file, "2 1 0 5").unwrap();
        for tag in tags {
            writeln!(file, "{tag}").unwrap();
        }
        writeln!(file, "{stray}").unwrap();
        writeln!(file, "0 0 0").unwrap();
        writeln!(file, "1 0 0").unwrap();
        writeln!(file, "1 1 0").unwrap();
        writeln!(file, "0 1 0").unwrap();
        writeln!(file, "0.5 2 0").unwrap();
        writeln!(file, "$EndNodes").unwrap();
        writeln!(file, "$Elements").unwrap();
        writeln!(file, "6 7 1 7").unwrap();
        // The four boundary segments of Γ, one curve block each.
        writeln!(file, "1 1 1 1").unwrap();
        writeln!(file, "1 {n1} {n2}").unwrap();
        writeln!(file, "1 2 1 1").unwrap();
        writeln!(file, "2 {n2} {n3}").unwrap();
        writeln!(file, "1 3 1 1").unwrap();
        writeln!(file, "3 {n3} {n4}").unwrap();
        writeln!(file, "1 4 1 1").unwrap();
        writeln!(file, "4 {n4} {n1}").unwrap();
        // The stray segment (Γbad).
        writeln!(file, "1 5 1 1").unwrap();
        writeln!(file, "5 {n1} {stray}").unwrap();
        // The two triangles.
        writeln!(file, "2 1 2 2").unwrap();
        writeln!(file, "6 {n1} {n2} {n3}").unwrap();
        writeln!(file, "7 {n3} {n4} {n1}").unwrap();
        writeln!(file, "$EndElements").unwrap();
        file
    }

    fn write_two_triangle_mesh() -> tempfile::NamedTempFile {
        write_two_triangle_mesh_with_tags([1, 2, 3, 4])
    }

    #[test]
    fn crtri3_boundary_reads_trace_elements() {
        // The four outer edges of the square, each traced from its triangle:
        // triangles are (0,1,2) and (2,3,0), so the boundary edges are
        //   (0,1) → triangle 0, facet 2   (2,3) → triangle 1, facet 2
        //   (1,2) → triangle 0, facet 0   (0,3) → triangle 1, facet 0
        let file = write_two_triangle_mesh();
        let mesh = read_gmsh_with::<2>(
            file.path(),
            GmshOptions {
                renumber: true,
                edges: true,
                faces: false,
            },
        )
        .unwrap();

        let crs = mesh.crtri3_elements("Ω").unwrap();
        let facets = mesh.crtri3_boundary("Γ").unwrap();
        assert_eq!(facets.len(), 4);

        let expected = [
            (0_usize, 2_usize), // edge (0,1) of parent 0
            (0, 0),             // edge (1,2) of parent 0
            (1, 2),             // edge (2,3) of parent 1
            (1, 0),             // edge (0,3) of parent 1
        ];
        for (facet, &(parent, local)) in facets.iter().zip(&expected) {
            assert_eq!(facet.local_facet(), local);
            assert_eq!(
                facet.id(),
                crs[parent].id(),
                "the trace dofs are the parent's facet dofs"
            );
            // jacobe = |facet| / 2 = 1 / 2 for unit edges.
            assert!((facet.jacobe([0.0]) - 0.5).abs() < 1e-12);
            // The facet's own dof is 1 on the whole facet.
            for ξ in [-1.0, 0.0, 1.0] {
                assert!((facet.shape([ξ])[local] - 1.0).abs() < 1e-12);
            }
        }

        // The trace element is usable with the boundary operators: a penalty on a
        // constant datum loads only the facet's own dof.
        let (triplets, doublets) =
            crate::operation::heat::temperature_penalty(&facets, 10.0, |_| 2.0);
        let mut load = [0.0_f64; 6];
        for d in &doublets {
            load[d.index] += d.value;
        }
        for facet in &facets {
            let local = facet.id()[facet.local_facet()];
            // |e| = 1, datum 2 → α |e| 2 = 20 on the facet's own dof.
            assert!(
                (load[local] - 20.0).abs() < 1e-12,
                "dof {local}: {}",
                load[local]
            );
        }
        assert_eq!(triplets.len() % 4, 0); // three basis functions × two segments? sanity only
    }

    #[test]
    fn crtri3_patchtest_with_dirichlet_penalty() {
        // End-to-end: CR Poisson on the patch-test square with a linear exact
        // solution u = 1 + 2x − 3y, whose CR interpolant is u itself (generalized
        // patch test). Dirichlet data are imposed on the trace elements with a
        // penalty, so the discrete solution must reproduce the edge averages.
        use crate::operation::heat::{conduction_stiffness, temperature_penalty};
        use crate::sparse_vector::try_new_from_doublet;
        use faer::Side;
        use faer::prelude::*;
        use faer::sparse::SparseColMat;

        let mesh = read_gmsh_with::<2>(
            "examples/msh/patchtest_tri3.msh",
            GmshOptions {
                renumber: true,
                edges: true,
                faces: false,
            },
        )
        .unwrap();
        let tris = mesh.crtri3_elements("Ω").unwrap();
        let boundary = mesh.crtri3_boundary("Γᵍ").unwrap();
        let n_dofs = mesh.edge_numbering().unwrap().len();

        let exact = |x: [f64; 2]| 1.0 + 2.0 * x[0] - 3.0 * x[1];
        let penalty = 1e7;

        let mut triplets = conduction_stiffness(&tris, 1.0);
        let (penalty_k, penalty_f) = temperature_penalty(&boundary, penalty, exact);
        triplets.extend(penalty_k);

        let k =
            SparseColMat::<usize, f64>::try_new_from_triplets(n_dofs, n_dofs, &triplets).unwrap();
        let f = try_new_from_doublet(n_dofs, &penalty_f).unwrap();
        let d = k.sp_cholesky(Side::Lower).unwrap().solve(&f);

        // Each edge dof is the edge average of a linear field, i.e. its midpoint value.
        let mut worst = 0.0_f64;
        for cr in &tris {
            for (a, id) in cr.id().iter().enumerate() {
                let midpoint = cr.nodes()[a].coordinates;
                worst = worst.max((d[*id] - exact(midpoint)).abs());
            }
        }
        assert!(worst < 1e-6, "patch test error {worst:e}");
    }

    /// Writes a 3D mesh of two tetrahedra sharing the face `(2, 3, 4)`:
    /// `t1 = (1, 2, 3, 4)` and `t2 = (5, 2, 3, 4)` with node 5 at `(-1,-1,-1)`.
    ///
    /// Groups: `Ω` (the two tetrahedra, physical tag 1) and `Γ` (the six outer
    /// triangles, physical tag 2).
    fn write_two_tet_mesh() -> tempfile::NamedTempFile {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "$MeshFormat").unwrap();
        writeln!(file, "4.1 0 8").unwrap();
        writeln!(file, "$EndMeshFormat").unwrap();
        writeln!(file, "$PhysicalNames").unwrap();
        writeln!(file, "2").unwrap();
        writeln!(file, "3 1 \"Ω\"").unwrap();
        writeln!(file, "2 2 \"Γ\"").unwrap();
        writeln!(file, "$EndPhysicalNames").unwrap();
        writeln!(file, "$Entities").unwrap();
        writeln!(file, "0 0 6 2").unwrap();
        // Six boundary triangles (physical tag 2), one surface entity each.
        for surface in 1..=6 {
            writeln!(file, "{surface} 0 0 0 1 1 1 1 2 0").unwrap();
        }
        // Two volume entities (physical tag 1).
        writeln!(file, "1 0 0 0 1 1 1 1 1 0").unwrap();
        writeln!(file, "2 0 0 0 1 1 1 1 1 0").unwrap();
        writeln!(file, "$EndEntities").unwrap();
        writeln!(file, "$Nodes").unwrap();
        writeln!(file, "1 5 1 5").unwrap();
        writeln!(file, "3 1 0 5").unwrap();
        for tag in 1..=5 {
            writeln!(file, "{tag}").unwrap();
        }
        writeln!(file, "0 0 0").unwrap();
        writeln!(file, "1 0 0").unwrap();
        writeln!(file, "0 1 0").unwrap();
        writeln!(file, "0 0 1").unwrap();
        writeln!(file, "-1 -1 -1").unwrap();
        writeln!(file, "$EndNodes").unwrap();
        writeln!(file, "$Elements").unwrap();
        writeln!(file, "9 9 1 9").unwrap();
        // Boundary triangles of Γ: three per tetrahedron, the shared face excluded.
        let boundary = [
            [2, 3, 4],
            [1, 3, 4],
            [1, 2, 4],
            [1, 2, 3],
            [5, 3, 4],
            [5, 2, 4],
            [5, 2, 3],
        ];
        for (index, nodes) in boundary.iter().enumerate().take(6) {
            writeln!(file, "2 {} 2 1", index + 1).unwrap();
            writeln!(file, "{} {} {} {}", index + 1, nodes[0], nodes[1], nodes[2]).unwrap();
        }
        writeln!(file, "2 7 2 1").unwrap();
        writeln!(
            file,
            "7 {} {} {}",
            boundary[6][0], boundary[6][1], boundary[6][2]
        )
        .unwrap();
        // The two tetrahedra.
        writeln!(file, "3 1 4 1").unwrap();
        writeln!(file, "8 1 2 3 4").unwrap();
        writeln!(file, "3 2 4 1").unwrap();
        writeln!(file, "8 5 2 3 4").unwrap();
        writeln!(file, "$EndElements").unwrap();
        file
    }

    #[test]
    fn crtet4_elements_and_boundary_share_face_dofs() {
        let file = write_two_tet_mesh();
        let mesh = read_gmsh_with::<3>(
            file.path(),
            GmshOptions {
                renumber: true,
                edges: false,
                faces: true,
            },
        )
        .unwrap();

        let tets = mesh.crtet4_elements("Ω").unwrap();
        assert_eq!(tets.len(), 2);
        let faces = mesh.face_numbering().unwrap();
        // Euler for a ball of two tets glued on one face: F = 4 + 4 − 1 = 7.
        assert_eq!(faces.len(), 7);
        assert_eq!(faces.offset(), 0);

        // All face dofs are distinct and inside `0..7`; the shared face (node ids
        // 1, 2, 3, i.e. GMSH tags 2, 3, 4) is one dof used by both tetrahedra.
        let mut ids: Vec<usize> = tets.iter().flat_map(|tet| tet.id()).collect();
        ids.sort_unstable();
        let all = ids.clone();
        ids.dedup();
        assert_eq!(ids.len(), 7);
        assert!(ids.iter().all(|&id| id < 7));
        let shared = faces.get(1, 2, 3).unwrap();
        assert!(tets[0].id().contains(&shared));
        assert!(tets[1].id().contains(&shared));
        assert_eq!(all.iter().filter(|&&id| id == shared).count(), 2);

        // The trace elements carry the parent's four face dofs.
        let boundary = mesh.crtet4_boundary("Γ").unwrap();
        assert_eq!(boundary.len(), 6);
        for facet in &boundary {
            assert!(facet.id().iter().all(|id| *id < 7));
            assert!(facet.local_facet() < 4);
        }
        // The trace of a boundary triangle belongs to exactly one tetrahedron.
        for facet in &boundary {
            let parents: Vec<&CRTet4<3>> =
                tets.iter().filter(|tet| tet.id() == facet.id()).collect();
            assert_eq!(
                parents.len(),
                1,
                "a boundary trace must match one tetrahedron"
            );
            // jacobe = 2 × triangle area, and the facet's own dof is 1 everywhere.
            assert!(facet.jacobe([1.0 / 3.0, 1.0 / 3.0]) > 0.0);
            assert!((facet.shape([0.4, 0.3])[facet.local_facet()] - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn crtet4_elements_require_the_faces_option() {
        let file = write_two_tet_mesh();
        let mesh = read_gmsh_with::<3>(
            file.path(),
            GmshOptions {
                renumber: true,
                edges: false,
                faces: false,
            },
        )
        .unwrap();
        assert!(mesh.face_numbering().is_none());
        assert!(matches!(
            mesh.crtet4_elements("Ω"),
            Err(MeshError::FaceNumberingMissing)
        ));
        assert!(matches!(
            mesh.crtet4_boundary("Γ"),
            Err(MeshError::FaceNumberingMissing)
        ));
        // With `faces` enabled the wrong group types are rejected.
        let mesh = read_gmsh_with::<3>(
            file.path(),
            GmshOptions {
                renumber: true,
                edges: false,
                faces: true,
            },
        )
        .unwrap();
        assert!(matches!(
            mesh.crtet4_elements("Γ"),
            Err(MeshError::ElementTypeMismatch { .. })
        ));
        assert!(matches!(
            mesh.crtet4_boundary("Ω"),
            Err(MeshError::ElementTypeMismatch { .. })
        ));
    }

    #[test]
    fn tet4_elements_are_read_from_gmsh() {
        let file = write_two_tet_mesh();
        let mesh = read_gmsh::<3>(file.path()).unwrap();
        let tets: Vec<Tet4<3>> = mesh.elements("Ω").unwrap();
        assert_eq!(tets.len(), 2);
        assert_eq!(tets[0].id(), [0, 1, 2, 3]);
        assert_eq!(tets[1].id(), [4, 1, 2, 3]);
        assert!((tets[0].volume() - 1.0 / 6.0).abs() < 1e-12);
    }

    #[test]
    fn crtri3_boundary_rejects_non_edge_segments() {
        let file = write_two_triangle_mesh();
        let mesh = read_gmsh_with::<2>(
            file.path(),
            GmshOptions {
                renumber: true,
                edges: true,
                faces: false,
            },
        )
        .unwrap();
        assert!(matches!(
            mesh.crtri3_boundary("Γbad"),
            Err(MeshError::MissingEdge { .. })
        ));
        assert!(matches!(
            mesh.crtri3_boundary("Ω"),
            Err(MeshError::ElementTypeMismatch { .. })
        ));
    }

    #[test]
    fn crtri3_boundary_requires_the_edges_option() {
        let file = write_two_triangle_mesh();
        let mesh = read_gmsh_with::<2>(
            file.path(),
            GmshOptions {
                renumber: true,
                edges: false,
                faces: false,
            },
        )
        .unwrap();
        assert!(matches!(
            mesh.crtri3_boundary("Γ"),
            Err(MeshError::EdgeNumberingMissing)
        ));
    }

    #[test]
    fn edge_table_is_orientation_independent() {
        let mut edges = EdgeTable::new(7);
        assert!(edges.is_empty());
        assert_eq!(edges.offset(), 7);
        assert_eq!(edges.id(4, 2), 7);
        assert_eq!(edges.id(2, 4), 7, "the same edge either way round");
        assert_eq!(edges.get(4, 2), Some(7));
        assert_eq!(edges.id(2, 5), 8);
        assert_eq!(edges.get(0, 1), None);
        assert_eq!(edges.len(), 2);
    }

    #[test]
    fn crtri3_elements_share_interface_dofs() {
        // Two triangles sharing the diagonal: the shared edge must get one single
        // CR dof, so that both elements reference the same freedom.
        let file = write_two_triangle_mesh();
        let options = GmshOptions {
            renumber: true,
            edges: true,
            faces: false,
        };
        let mesh = read_gmsh_with::<2>(file.path(), options).unwrap();
        let edges = mesh.edge_numbering().unwrap();
        assert_eq!(edges.offset(), 0, "edge dofs are numbered from 0");

        let crs = mesh.crtri3_elements("Ω").unwrap();
        assert_eq!(crs.len(), 2);
        // Euler: V − E + F = 1 for a simply connected triangulation → E = 5.
        assert_eq!(edges.len(), 5);

        // Local edge order is (opposite vertex 1, 2, 3) = ((1,2), (0,2), (0,1)).
        // Triangle 1 = (0,1,2) and triangle 2 = (2,3,0) share the edge (0,2), which
        // is local edge 2 in both.
        assert_eq!(crs[0].id()[1], crs[1].id()[1]);
        assert_eq!(edges.get(0, 2), Some(crs[0].id()[1]));
        // All other edge dofs are distinct and lie in `0..5`.
        let mut ids: Vec<usize> = crs.iter().flat_map(|cr| cr.id()).collect();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), 5);
        assert!(ids.iter().all(|&id| (0..5).contains(&id)));

        // The midpoints are the edge midpoints of the physical mesh.
        assert_eq!(crs[0].nodes()[0].coordinates, [1.0, 0.5]); // edge (1,2)
        assert_eq!(crs[0].nodes()[1].coordinates, [0.5, 0.5]); // edge (0,2)
        assert_eq!(crs[0].nodes()[2].coordinates, [0.5, 0.0]); // edge (0,1)
        assert_eq!(
            crs[0].triangle_vertices(),
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
        );
    }

    #[test]
    fn crtri3_elements_require_the_edges_option() {
        // Without `GmshOptions::edges` the mesh carries no edge numbering, and the
        // explicit edge table path is the way in.
        let file = write_two_triangle_mesh();
        let mesh = read_gmsh_with::<2>(
            file.path(),
            GmshOptions {
                renumber: true,
                edges: false,
                faces: false,
            },
        )
        .unwrap();
        assert!(mesh.edge_numbering().is_none());
        assert!(matches!(
            mesh.crtri3_elements("Ω"),
            Err(MeshError::EdgeNumberingMissing)
        ));

        // An explicit table gives the same ids (same visiting order).
        let table = mesh.node_table().unwrap();
        let mut edges = mesh.edge_table();
        assert_eq!(edges.offset(), 0, "the edge numbering starts at 0");
        let crs = mesh
            .crtri3_elements_with_table("Ω", &table, &mut edges)
            .unwrap();
        assert_eq!(crs.len(), 2);
        assert_eq!(edges.len(), 5);

        let read_time = read_gmsh_with::<2>(
            file.path(),
            GmshOptions {
                renumber: true,
                edges: true,
                faces: false,
            },
        )
        .unwrap();
        assert_eq!(read_time.crtri3_elements("Ω").unwrap()[0].id(), crs[0].id());
    }

    #[test]
    fn crtri3_elements_from_patchtest_mesh() {
        // A real mesh: 288 triangles over a square. Every triangle gets three CR
        // dofs, all edges are shared exactly once, and the edge count matches the
        // number of distinct vertex pairs.
        let mesh = read_gmsh_with::<2>(
            "examples/msh/patchtest_tri3.msh",
            GmshOptions {
                renumber: true,
                edges: false,
                faces: false,
            },
        )
        .unwrap();
        let table = mesh.node_table().unwrap();
        let tris: Vec<Tri3<2>> = mesh.elements("Ω").unwrap();

        let mut edges = mesh.edge_table();
        let crs = mesh
            .crtri3_elements_with_table("Ω", &table, &mut edges)
            .unwrap();
        assert_eq!(crs.len(), tris.len());

        let mut distinct: std::collections::HashSet<(usize, usize)> =
            std::collections::HashSet::new();
        for tri in &tris {
            let [a, b, c] = tri.id();
            for (i, j) in [(a, b), (b, c), (a, c)] {
                distinct.insert((i.min(j), i.max(j)));
            }
        }
        assert_eq!(edges.len(), distinct.len());
        assert_eq!(edges.len(), table.len() + tris.len() - 1, "Euler formula");

        for cr in &crs {
            for id in cr.id() {
                assert!(
                    id < edges.len(),
                    "edge dof {id} outside the edge range 0..{}",
                    edges.len()
                );
            }
        }
        // Reading the same group again reuses the existing ids.
        let before = edges.len();
        let again = mesh
            .crtri3_elements_by_tag_with_table(
                mesh.physical_group_by_name("Ω").unwrap().tag,
                &table,
                &mut edges,
            )
            .unwrap();
        assert_eq!(again.len(), crs.len());
        assert_eq!(again[0].id(), crs[0].id());
        assert_eq!(edges.len(), before);

        // With `edges: true` the numbering is created while reading, spans all
        // triangle groups and requires no table argument at all.
        let prepared = read_gmsh_with::<2>(
            "examples/msh/patchtest_tri3.msh",
            GmshOptions {
                renumber: true,
                edges: true,
                faces: false,
            },
        )
        .unwrap();
        let numbering = prepared.edge_numbering().unwrap();
        assert_eq!(numbering.offset(), 0);
        assert_eq!(numbering.len(), edges.len());
        assert!(!numbering.is_empty());
        let direct = prepared.crtri3_elements("Ω").unwrap();
        assert_eq!(direct.len(), tris.len());
        for (a, b) in direct.iter().zip(&crs) {
            assert_eq!(a.id(), b.id());
        }
    }

    #[test]
    fn crtri3_elements_reject_non_triangle_groups() {
        let mesh = read_gmsh_with::<2>(
            "examples/msh/patchtest_tri3.msh",
            GmshOptions {
                renumber: true,
                edges: true,
                faces: false,
            },
        )
        .unwrap();

        assert!(matches!(
            mesh.crtri3_elements("Γᵍ"),
            Err(MeshError::ElementTypeMismatch { .. })
        ));
    }

    #[test]
    fn crtri3_elements_accept_a_block_offset() {
        // The edge numbering may be placed after another field's dofs (mixed
        // systems); the read-time numbering keeps `0..n_edges`, so both layouts are
        // available and differ by a constant shift only.
        let file = write_two_triangle_mesh();
        let mesh = read_gmsh_with::<2>(
            file.path(),
            GmshOptions {
                renumber: true,
                edges: false,
                faces: false,
            },
        )
        .unwrap();
        let table = mesh.node_table().unwrap();

        let mut shifted = EdgeTable::new(100);
        let crs = mesh
            .crtri3_elements_with_table("Ω", &table, &mut shifted)
            .unwrap();
        let ids: Vec<usize> = crs.iter().flat_map(|cr| cr.id()).collect();
        assert!(ids.iter().all(|&id| (100..105).contains(&id)));
        assert_eq!(shifted.offset(), 100);
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
