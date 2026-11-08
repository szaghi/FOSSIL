do i = 1, mesh%get_patches_number()
   n = count(mesh%facet_patch(facets()) == i)
   if (n > 0) print '(3A,I0,A)', 'patch ', mesh%patch_name(i), ': ', n, ' facets'   ! empty ones skipped
enddo
