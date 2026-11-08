surface = dragon
call surface%resize(factor=factor, respect_centroid=.true.)  ! about the centroid of the last analyze
call surface%analyze
