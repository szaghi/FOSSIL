import { withMermaid } from 'vitepress-plugin-mermaid'
import apiSidebar from '../api/_sidebar.json'

// one sidebar for every page but the API, in reading order: the "previous" and "next" links at the bottom of a page
// follow it, so the documentation reads from the first page to the last
const docs = [
  {
    text: 'Start here',
    items: [
      { text: 'Introduction', link: '/guide/' },
      { text: 'Installation', link: '/guide/installation' },
    ],
  },
  {
    text: 'Tutorial',
    items: [
      { text: 'Overview', link: '/manual/' },
      { text: '1. First surface', link: '/manual/tutorial/01-first-surface' },
      { text: '2. Repair', link: '/manual/tutorial/02-repair' },
      { text: '3. Manipulation', link: '/manual/tutorial/03-manipulation' },
      { text: '4. Distance and sign', link: '/manual/tutorial/04-distance' },
      { text: '5. Batched queries', link: '/manual/tutorial/05-batched' },
      { text: '6. Distance fields', link: '/manual/tutorial/06-distance-field' },
      { text: '7. 2-D curves', link: '/manual/tutorial/07-curves' },
      { text: '8. Immersed-boundary queries', link: '/manual/tutorial/08-ib-queries' },
      { text: '9. AMR tagging', link: '/manual/tutorial/09-amr-tagging' },
      { text: '10. Cut cells', link: '/manual/tutorial/10-cut-cells' },
      { text: '11. Moving bodies', link: '/manual/tutorial/11-moving-bodies' },
      { text: '12. On the GPU', link: '/manual/tutorial/12-gpu' },
    ],
  },
  {
    text: 'Recipes',
    items: [{ text: 'Cookbook', link: '/manual/cookbook' }],
  },
  {
    text: 'Reference: surfaces',
    collapsed: false,
    items: [
      { text: 'Feature map', link: '/guide/features' },
      { text: 'I/O, analysis, manipulation', link: '/guide/usage' },
      { text: 'Hole filling', link: '/guide/advanced/hole-filling' },
      { text: 'Self-intersection', link: '/guide/advanced/self-intersection' },
      { text: 'Boolean operations', link: '/guide/advanced/booleans' },
      { text: 'Mesh decimation', link: '/guide/advanced/decimate' },
      { text: 'Isotropic remeshing', link: '/guide/advanced/isotropic-remesh' },
      { text: 'Mesh smoothing', link: '/guide/advanced/smoothing' },
      { text: 'Alpha wrapping', link: '/guide/advanced/alpha-wrap' },
      { text: 'Marching cubes', link: '/guide/advanced/marching-cubes' },
      { text: 'SDF segmentation', link: '/guide/advanced/sdf-segmentation' },
      { text: 'Cotangent Laplacian', link: '/guide/advanced/cotangent-laplacian' },
      { text: 'Per-vertex curvature', link: '/guide/advanced/curvature' },
      { text: 'Mass properties', link: '/guide/advanced/mass-properties' },
      { text: 'Surface deviation', link: '/guide/advanced/surface-deviation' },
      { text: 'Feature edges', link: '/guide/advanced/feature-edges' },
      { text: 'Thickness and gaps', link: '/guide/advanced/thickness' },
      { text: 'Named patches', link: '/guide/advanced/patches' },
      { text: 'Body–body distance', link: '/guide/advanced/body-distance' },
    ],
  },
  {
    text: 'Reference: distance and signs',
    collapsed: false,
    items: [
      { text: 'Signed-distance fields', link: '/guide/advanced/distance-field' },
      { text: 'Generalized winding number', link: '/guide/advanced/winding-number' },
      { text: 'Ray queries', link: '/guide/advanced/ray-queries' },
      { text: 'Sign algorithms benchmark', link: '/guide/advanced/sign-benchmark' },
    ],
  },
  {
    text: 'Reference: immersed boundaries',
    collapsed: false,
    items: [
      { text: '2-D geometry (curves)', link: '/guide/advanced/2d-geometry' },
      { text: 'Immersed-boundary queries', link: '/guide/advanced/ib-queries' },
      { text: 'AMR tagging', link: '/guide/advanced/amr-tagging' },
      { text: 'Cut-cell geometry', link: '/guide/advanced/cut-cells' },
      { text: 'Cut cells: solver support', link: '/guide/advanced/cut-solver' },
      { text: 'Moving bodies', link: '/guide/advanced/moving-bodies' },
    ],
  },
  {
    text: 'Reference: GPU',
    collapsed: false,
    items: [{ text: 'GPU queries', link: '/guide/advanced/device-distance' }],
  },
  {
    text: 'Reference: API companion',
    collapsed: true,
    items: [
      { text: 'Overview', link: '/guide/api/' },
      { text: 'Constants', link: '/guide/api/constants' },
      { text: 'surface_stl_object', link: '/guide/api/surface-stl-object' },
      { text: 'facet_object', link: '/guide/api/facet-object' },
      { text: 'aabb_tree_object', link: '/guide/api/aabb-tree-object' },
      { text: 'curve_object', link: '/guide/api/curve-object' },
    ],
  },
  {
    text: 'Project',
    items: [
      { text: 'Changelog', link: '/guide/changelog' },
      { text: 'Contributing', link: '/guide/contributing' },
    ],
  },
]

export default withMermaid({
  title: 'FOSSIL',
  description: 'FOrtran Surfaces, Solids, Immersed boundaries and Level sets: the geometry of complex bodies for simulation codes, in pure Fortran.',
  base: '/FOSSIL/',
  head: [['link', { rel: 'icon', type: 'image/svg+xml', href: '/FOSSIL/favicon.svg' }]],
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
    logo: { light: '/logo.svg', dark: '/logo-dark.svg', alt: 'FOSSIL logo' },
    socialLinks: [{ icon: 'github', link: 'https://github.com/szaghi/FOSSIL' }],
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Start here', link: '/guide/', activeMatch: '^/guide/(index|installation)?$' },
      { text: 'Tutorial', link: '/manual/tutorial/01-first-surface', activeMatch: '^/manual/(index|tutorial/)' },
      { text: 'Cookbook', link: '/manual/cookbook', activeMatch: '^/manual/cookbook' },
      { text: 'Reference', link: '/guide/features', activeMatch: '^/guide/(features|usage|advanced/|api/)' },
      { text: 'API', link: '/api/', activeMatch: '^/api/' },
      {
        text: 'Project',
        items: [
          { text: 'Changelog', link: '/guide/changelog' },
          { text: 'Contributing', link: '/guide/contributing' },
        ],
      },
    ],
    sidebar: {
      '/guide/': docs,
      '/manual/': docs,
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
    build: {
      target: 'es2022',
    },
    optimizeDeps: {
      include: ['mermaid'],
    },
  },
})
