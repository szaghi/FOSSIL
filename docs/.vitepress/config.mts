import { withMermaid } from 'vitepress-plugin-mermaid'
import apiSidebar from '../api/_sidebar.json'

export default withMermaid({
  title: 'FOSSIL Documentation',
  base: '/FOSSIL/',
  markdown: {
    math: true,
    languages: ['fortran-free-form', 'fortran-fixed-form'],
    languageAlias: {
      'fortran': 'fortran-free-form',
      'f90': 'fortran-free-form',
      'f95': 'fortran-free-form',
      'f03': 'fortran-free-form',
      'f08': 'fortran-free-form',
      'f77': 'fortran-fixed-form',
    },
  },
  themeConfig: {
    nav: [
      { text: 'Home', link: '/' },
      {
        text: 'Guide',
        items: [
          { text: 'About',             link: '/guide/' },
          { text: 'Features',          link: '/guide/features' },
          { text: 'Installation',      link: '/guide/installation' },
          { text: 'Usage',             link: '/guide/usage' },
          { text: 'API Companion Guide', link: '/guide/api/' },
          { text: 'Contributing',      link: '/guide/contributing' },
          { text: 'Coverage Analysis', link: '/guide/coverage-analysis' },
          { text: 'Changelog',         link: '/guide/changelog' },
        ],
      },
      { text: 'API', link: '/api/' },
      { text: 'GitHub', link: 'https://github.com/szaghi/FOSSIL' },
    ],
    sidebar: {
      '/guide/': [
        {
          text: 'Introduction',
          items: [
            { text: 'About',    link: '/guide/' },
            { text: 'Features', link: '/guide/features' },
          ],
        },
        {
          text: 'Getting Started',
          items: [
            { text: 'Installation',  link: '/guide/installation' },
            { text: 'Usage',         link: '/guide/usage' },
          ],
        },
        {
          text: 'API Companion Guide',
          items: [
            { text: 'Overview',             link: '/guide/api/' },
            { text: 'Constants',            link: '/guide/api/constants' },
            { text: 'surface_stl_object',   link: '/guide/api/surface-stl-object' },
            { text: 'facet_object',         link: '/guide/api/facet-object' },
            { text: 'aabb_tree_object',     link: '/guide/api/aabb-tree-object' },
          ],
        },
        {
          text: 'Advanced Features',
          items: [
            { text: 'Overview',                  link: '/guide/advanced/' },
            { text: 'Boolean operations',        link: '/guide/advanced/booleans' },
            { text: 'Self-intersection',         link: '/guide/advanced/self-intersection' },
            { text: 'Mesh decimation',           link: '/guide/advanced/decimate' },
            { text: 'Generalized winding number', link: '/guide/advanced/winding-number' },
            { text: 'Marching cubes',            link: '/guide/advanced/marching-cubes' },
            { text: 'Alpha wrapping',            link: '/guide/advanced/alpha-wrap' },
            { text: 'Isotropic remeshing',       link: '/guide/advanced/isotropic-remesh' },
            { text: 'SDF segmentation',          link: '/guide/advanced/sdf-segmentation' },
            { text: 'Ray queries',               link: '/guide/advanced/ray-queries' },
          ],
        },
        {
          text: 'Project',
          items: [
            { text: 'Contributing',      link: '/guide/contributing' },
            { text: 'Coverage Analysis', link: '/guide/coverage-analysis' },
            { text: 'Changelog',         link: '/guide/changelog' },
          ],
        },
      ],
      '/api/': [
        {
          text: 'API Reference',
          items: [
            { text: 'Overview', link: '/api/' },
          ],
        },
        ...apiSidebar,
      ],
    },
    search: {
      provider: 'local',
    },
  },
  mermaid: {},
  vite: {
    optimizeDeps: {
      include: ['mermaid'],
    },
  },
})
