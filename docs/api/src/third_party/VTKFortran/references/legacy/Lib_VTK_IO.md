---
title: Lib_VTK_IO
---

# Lib_VTK_IO

> Pure Fortran (2003+) library to write and read data conforming the VTK standard
{!README-Lib_VTK_IO.md!}

**Source**: `src/third_party/VTKFortran/references/legacy/Lib_VTK_IO.f90`

**Dependencies**

```mermaid
graph LR
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_CON["Lib_VTK_IO_CON"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_CON_XML["Lib_VTK_IO_CON_XML"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_DAT_VAR["Lib_VTK_IO_DAT_VAR"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_DAT_VAR_XML["Lib_VTK_IO_DAT_VAR_XML"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_END["Lib_VTK_IO_END"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_END_XML["Lib_VTK_IO_END_XML"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_FLD_XML["Lib_VTK_IO_FLD_XML"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_GEO["Lib_VTK_IO_GEO"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_GEO_XML["Lib_VTK_IO_GEO_XML"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_INI["Lib_VTK_IO_INI"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_INI_XML["Lib_VTK_IO_INI_XML"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_PVD_XML["Lib_VTK_IO_PVD_XML"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_PVTK_XML["Lib_VTK_IO_PVTK_XML"]
  Lib_VTK_IO["Lib_VTK_IO"] --> Lib_VTK_IO_VTM_XML["Lib_VTK_IO_VTM_XML"]
```
