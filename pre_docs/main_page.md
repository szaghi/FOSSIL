project: FOSSIL
src_dir: ../src/lib
src_dir: ../src/tests
exclude_dir: ../src/third_party/PENF/src/lib
             ../src/third_party/PENF/src/tests
output_dir: ../docs/
project_github: https://github.com/szaghi/FOSSIL
summary: FOSSIL, FOrtran Stereo (si) Litography parser
author: Stefano Zaghi
github: https://github.com/szaghi
email: stefano.zaghi@gmail.com
md_extensions: markdown.extensions.toc(anchorlink=True)
               markdown.extensions.smarty(smart_quotes=False)
               markdown.extensions.extra
               markdown_checklist.extension
docmark: <
display: public
         protected
         private
source: true
warn: true
graph: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

{!README-FOSSIL.md!}
