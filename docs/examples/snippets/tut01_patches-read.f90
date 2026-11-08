call back%load_from_file(file_name='dragon-patches.stl', guess_format=.true.)
print '(A,I0)', 'patches read back ', back%get_patches_number()
do p = 1, back%get_patches_number()
   print '(A,I0,A,A8,I6,A)', 'patch ', p, '  ', back%patch_name(p), &
                             count(back%facet_patch(ids) == p), ' facets'
enddo
