print '(A,I0,A,6F9.5)', '  moments up to order 2 (', mono_count(2_I4P, 2_I4P), '), / cell area: ', &
                        cut%piece_moments(:, p) / 0.01_R8P
