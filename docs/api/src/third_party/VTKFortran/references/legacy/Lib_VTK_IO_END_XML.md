---
title: Lib_VTK_IO_END_XML
---

# Lib_VTK_IO_END_XML

> END_XML interface definition for Lib_VTK_IO.

**Source**: `src/third_party/VTKFortran/references/legacy/Lib_VTK_IO_END_XML.f90`

**Dependencies**

```mermaid
graph LR
  Lib_VTK_IO_END_XML["Lib_VTK_IO_END_XML"] --> Lib_VTK_IO_Back_End["Lib_VTK_IO_Back_End"]
  Lib_VTK_IO_END_XML["Lib_VTK_IO_END_XML"] --> befor64["befor64"]
  Lib_VTK_IO_END_XML["Lib_VTK_IO_END_XML"] --> penf["penf"]
```

## Contents

- [VTK_END_XML_WRITE](#vtk-end-xml-write)
- [VTK_END_XML_READ](#vtk-end-xml-read)

## Functions

### VTK_END_XML_WRITE

Function for finalizing the VTK-XML file.

### Usage
```fortran
 E_IO = VTK_END_XML()
```

**Returns**: integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables))

```fortran
function VTK_END_XML_WRITE(cf) result(E_IO)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `cf` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | inout | optional | Current file index (for concurrent files IO). |

**Call graph**

```mermaid
flowchart TD
  test_mpi["test_mpi"] --> VTK_END_XML_WRITE["VTK_END_XML_WRITE"]
  test_openmp["test_openmp"] --> VTK_END_XML_WRITE["VTK_END_XML_WRITE"]
  test_pstrg["test_pstrg"] --> VTK_END_XML_WRITE["VTK_END_XML_WRITE"]
  test_punst["test_punst"] --> VTK_END_XML_WRITE["VTK_END_XML_WRITE"]
  test_rect["test_rect"] --> VTK_END_XML_WRITE["VTK_END_XML_WRITE"]
  test_strg["test_strg"] --> VTK_END_XML_WRITE["VTK_END_XML_WRITE"]
  test_unst["test_unst"] --> VTK_END_XML_WRITE["VTK_END_XML_WRITE"]
  test_vtm["test_vtm"] --> VTK_END_XML_WRITE["VTK_END_XML_WRITE"]
  VTK_END_XML_WRITE["VTK_END_XML_WRITE"] --> b64_encode["b64_encode"]
  VTK_END_XML_WRITE["VTK_END_XML_WRITE"] --> pack_data["pack_data"]
  VTK_END_XML_WRITE["VTK_END_XML_WRITE"] --> str["str"]
  VTK_END_XML_WRITE["VTK_END_XML_WRITE"] --> vtk_update["vtk_update"]
  style VTK_END_XML_WRITE fill:#3e63dd,stroke:#99b,stroke-width:2px
```

### VTK_END_XML_READ

Function for close an opened VTK file.

**Returns**: integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables))

```fortran
function VTK_END_XML_READ(cf) result(E_IO)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `cf` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) |  | optional |  |
