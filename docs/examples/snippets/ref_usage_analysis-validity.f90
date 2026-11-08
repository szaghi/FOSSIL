if (.not. surface%is_watertight()) print '(A)', 'surface has boundary or non-manifold edges'
if (.not. surface%is_manifold())   print '(A)', 'surface has non-manifold edges (3+ incidences)'
if (.not. surface%is_volume()) then
   print '(A)', 'volume is not physically meaningful: do not trust get_volume()/get_centroid()'
endif
