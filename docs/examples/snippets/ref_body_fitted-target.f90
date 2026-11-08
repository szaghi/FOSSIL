call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.) ! chord x in [-0.25, 0.75], span z in [0, 1]
allocate(flap(0))
do i = 1, wing%get_facets_number() ! the rear facets become the patch 'flap'
   f => wing%facet_at(i)
   if (min(f%vertex(1)%x, f%vertex(2)%x, f%vertex(3)%x) >= 0.45_R8P) flap = [flap, i]
enddo
call wing%set_patch(flap, 'flap')
call wing%projection_target(target) ! sharp edges above 30 degrees, patch boundaries, corners
print '(A,I0,A,I0,A,I0)', 'target: patches ', target%np, ', constraint curves ', target%nc, ', corners ', target%nk
