---
title: Lib_VTK_IO_PVTK_XML
---

# Lib_VTK_IO_PVTK_XML

> PVTK_XML interface definitions for Lib_VTK_IO.

**Source**: `src/third_party/VTKFortran/references/legacy/Lib_VTK_IO_PVTK_XML.f90`

**Dependencies**

```mermaid
graph LR
  Lib_VTK_IO_PVTK_XML["Lib_VTK_IO_PVTK_XML"] --> Lib_VTK_IO_Back_End["Lib_VTK_IO_Back_End"]
  Lib_VTK_IO_PVTK_XML["Lib_VTK_IO_PVTK_XML"] --> befor64["befor64"]
  Lib_VTK_IO_PVTK_XML["Lib_VTK_IO_PVTK_XML"] --> penf["penf"]
```

## Contents

- [PVTK_INI_XML](#pvtk-ini-xml)
- [PVTK_GEO_XML](#pvtk-geo-xml)
- [PVTK_DAT_XML](#pvtk-dat-xml)
- [PVTK_VAR_XML](#pvtk-var-xml)
- [PVTK_END_XML](#pvtk-end-xml)

## Functions

### PVTK_INI_XML

Function for initializing parallel (partitioned) VTK-XML file.

**Returns**: integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables))

```fortran
function PVTK_INI_XML(filename, mesh_topology, tp, cf, nx1, nx2, ny1, ny2, nz1, nz2) result(E_IO)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `filename` | character(len=*) | in |  | File name. |
| `mesh_topology` | character(len=*) | in |  | Mesh topology. |
| `tp` | character(len=*) | in |  | Type of geometry representation (Float32, Float64, ecc). |
| `cf` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | out | optional | Current file index (for concurrent files IO). |
| `nx1` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Initial node of x axis. |
| `nx2` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Final node of x axis. |
| `ny1` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Initial node of y axis. |
| `ny2` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Final node of y axis. |
| `nz1` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Initial node of z axis. |
| `nz2` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Final node of z axis. |

**Call graph**

```mermaid
flowchart TD
  test_mpi["test_mpi"] --> PVTK_INI_XML["PVTK_INI_XML"]
  test_openmp["test_openmp"] --> PVTK_INI_XML["PVTK_INI_XML"]
  test_pstrg["test_pstrg"] --> PVTK_INI_XML["PVTK_INI_XML"]
  test_punst["test_punst"] --> PVTK_INI_XML["PVTK_INI_XML"]
  PVTK_INI_XML["PVTK_INI_XML"] --> Get_Unit["Get_Unit"]
  PVTK_INI_XML["PVTK_INI_XML"] --> b64_init["b64_init"]
  PVTK_INI_XML["PVTK_INI_XML"] --> penf_init["penf_init"]
  PVTK_INI_XML["PVTK_INI_XML"] --> str["str"]
  PVTK_INI_XML["PVTK_INI_XML"] --> vtk_update["vtk_update"]
  style PVTK_INI_XML fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### PVTK_GEO_XML

Function for saving piece geometry source for parallel (partitioned) VTK-XML file.

**Returns**: integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables))

```fortran
function PVTK_GEO_XML(source, cf, nx1, nx2, ny1, ny2, nz1, nz2) result(E_IO)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `source` | character(len=*) | in |  | Source file name containing the piece data. |
| `cf` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Current file index (for concurrent files IO). |
| `nx1` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Initial node of x axis. |
| `nx2` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Final node of x axis. |
| `ny1` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Initial node of y axis. |
| `ny2` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Final node of y axis. |
| `nz1` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Initial node of z axis. |
| `nz2` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Final node of z axis. |

**Call graph**

```mermaid
flowchart TD
  test_mpi["test_mpi"] --> PVTK_GEO_XML["PVTK_GEO_XML"]
  test_openmp["test_openmp"] --> PVTK_GEO_XML["PVTK_GEO_XML"]
  test_pstrg["test_pstrg"] --> PVTK_GEO_XML["PVTK_GEO_XML"]
  test_punst["test_punst"] --> PVTK_GEO_XML["PVTK_GEO_XML"]
  PVTK_GEO_XML["PVTK_GEO_XML"] --> str["str"]
  style PVTK_GEO_XML fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### PVTK_DAT_XML

Function for initializing/finalizing the saving of data associated to the mesh.

 Function that **must** be called before saving the data related to geometric mesh, this function initializes the
 saving of data variables indicating the *type* (node or cell centered) of variables that will be saved.

**Returns**: integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables))

```fortran
function PVTK_DAT_XML(var_location, var_block_action, cf) result(E_IO)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `var_location` | character(len=*) | in |  | Location of saving variables: CELL or NODE centered. |
| `var_block_action` | character(len=*) | in |  | Variables block action: OPEN or CLOSE block. |
| `cf` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Current file index (for concurrent files IO). |

**Call graph**

```mermaid
flowchart TD
  test_mpi["test_mpi"] --> PVTK_DAT_XML["PVTK_DAT_XML"]
  test_openmp["test_openmp"] --> PVTK_DAT_XML["PVTK_DAT_XML"]
  test_pstrg["test_pstrg"] --> PVTK_DAT_XML["PVTK_DAT_XML"]
  test_punst["test_punst"] --> PVTK_DAT_XML["PVTK_DAT_XML"]
  PVTK_DAT_XML["PVTK_DAT_XML"] --> Upper_Case["Upper_Case"]
  style PVTK_DAT_XML fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### PVTK_VAR_XML

Function for saving variable associated to nodes or cells geometry.

**Returns**: integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables))

```fortran
function PVTK_VAR_XML(varname, tp, cf, Nc) result(E_IO)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `varname` | character(len=*) | in |  | Variable name. |
| `tp` | character(len=*) | in |  | Type of data representation (Float32, Float64, ecc). |
| `cf` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Current file index (for concurrent files IO). |
| `Nc` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in | optional | Number of components of variable. |

**Call graph**

```mermaid
flowchart TD
  test_mpi["test_mpi"] --> PVTK_VAR_XML["PVTK_VAR_XML"]
  test_openmp["test_openmp"] --> PVTK_VAR_XML["PVTK_VAR_XML"]
  test_pstrg["test_pstrg"] --> PVTK_VAR_XML["PVTK_VAR_XML"]
  test_punst["test_punst"] --> PVTK_VAR_XML["PVTK_VAR_XML"]
  PVTK_VAR_XML["PVTK_VAR_XML"] --> str["str"]
  style PVTK_VAR_XML fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### PVTK_END_XML

Function for finalizing the parallel (partitioned) VTK-XML file.

**Returns**: integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables))

```fortran
function PVTK_END_XML(cf) result(E_IO)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `cf` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | optional | Current file index (for concurrent files IO). |

**Call graph**

```mermaid
flowchart TD
  test_mpi["test_mpi"] --> PVTK_END_XML["PVTK_END_XML"]
  test_openmp["test_openmp"] --> PVTK_END_XML["PVTK_END_XML"]
  test_pstrg["test_pstrg"] --> PVTK_END_XML["PVTK_END_XML"]
  test_punst["test_punst"] --> PVTK_END_XML["PVTK_END_XML"]
  PVTK_END_XML["PVTK_END_XML"] --> vtk_update["vtk_update"]
  style PVTK_END_XML fill:#3e63dd,stroke:#99b,stroke-width:2px
```
