# Changelog
## [v1.0.9](https://github.com/szaghi/FLAP/tree/v1.0.9) (2026-02-22)
[Full Changelog](https://github.com/szaghi/FLAP/compare/v1.0.8...v1.0.9)
### Miscellaneous
- Update BeFoR64 submodule ([`6361436`](https://github.com/szaghi/FLAP/commit/6361436485e1ca3bffbd9d8e1e648d4e2b7af900))

## [v1.0.8](https://github.com/szaghi/FLAP/tree/v1.0.8) (2026-02-22)
[Full Changelog](https://github.com/szaghi/FLAP/compare/v1.0.7...v1.0.8)
### CI/CD
- Deploy pages even if coverage analysis fails ([`55553e5`](https://github.com/szaghi/FLAP/commit/55553e50cf1e854594e882600988d344eab3abe8))

## [v1.0.7](https://github.com/szaghi/FLAP/tree/v1.0.7) (2026-02-21)
[Full Changelog](https://github.com/szaghi/FLAP/compare/v1.0.6...v1.0.7)
### Bug fixes
- Update third party submodules ([`1ae2da2`](https://github.com/szaghi/FLAP/commit/1ae2da20a61ce724549634fcf7a56791ce7577aa))

### Miscellaneous
- Correct GHA ci ([`3463adb`](https://github.com/szaghi/FLAP/commit/3463adb783259a5a900c84ef1351b7c1a6e960e5))

### New features
- Add full VitePress guide and modernise project tooling ([`dbf2628`](https://github.com/szaghi/FLAP/commit/dbf26286947f54071a21208fa8642c15def3d273))

## [v1.0.6](https://github.com/szaghi/FLAP/tree/v1.0.6) (2022-07-06)
[Full Changelog](https://github.com/szaghi/FLAP/compare/v1.0.5...v1.0.6)
### Miscellaneous
- Merge tag 'v1.0.5' into develop

Fix issue[#2](https://github.com/szaghi/FLAP/issues/2)

Fix translate method of surface object. ([`0b65a87`](https://github.com/szaghi/FLAP/commit/0b65a876fca0086db8bad390234eb005e64283e0))
- Analysis not completed ([`33537ce`](https://github.com/szaghi/FLAP/commit/33537ce0e5833e969fbe93b90710d02aa976b6d2))
- Commit before sanite status ([`2067636`](https://github.com/szaghi/FLAP/commit/206763658e1287f3bee1636e4d1237efccfc986c))
- Merge branch 'feature/analize_and_improve_speed' into develop ([`24606ee`](https://github.com/szaghi/FLAP/commit/24606eee7aaa6aa517200ca27d924231d588f0c6))
- Update submodules ([`3f121cb`](https://github.com/szaghi/FLAP/commit/3f121cbd7b958c84476653250259779b1d6d62d4))
- Clean ([`0612912`](https://github.com/szaghi/FLAP/commit/0612912474fc5f7daf927ff48224f0b2f48bd4b9))
- Update submodules ([`411f11e`](https://github.com/szaghi/FLAP/commit/411f11e2a7e96ea0824ec77f0c0b9903e9b44192))
- Clean ([`95549e5`](https://github.com/szaghi/FLAP/commit/95549e56634dd2b012c44003ba16eb24e73dcab4))
- Update submodules ([`29f6508`](https://github.com/szaghi/FLAP/commit/29f6508946880759105561b31b4b0ecc847c06b2))
- Improve test distance API ([`b1906b1`](https://github.com/szaghi/FLAP/commit/b1906b1410b65451c2ab66d108a694144142b437))
- Add Immersed Boundary field app generator

Short description

Add a new program for computing Immersed Boundary field for Xall ([`7c9cea3`](https://github.com/szaghi/FLAP/commit/7c9cea3de456565051a788aa13e5da738ea21437))
- Add new IB generator ([`d5c0732`](https://github.com/szaghi/FLAP/commit/d5c0732f5de6b91fa4ebdcdce86bea6145dbfc58))
- Algorithm checking: AABB failing

Performed AABB algorithm check: the speedup is poor and the consistency
is failing: the computed distance values are varing at refinemnt level
varing.

Must be totally re-implemented. ([`0920641`](https://github.com/szaghi/FLAP/commit/09206413e276073b76f76e9402903fda04a13c61))
- Refresh on new laptop ([`cb2dd3d`](https://github.com/szaghi/FLAP/commit/cb2dd3d7374f90fd1b9f4e612d97c003e35a05f9))
- Minor change to enable NVFortran compilation ([`cff6dff`](https://github.com/szaghi/FLAP/commit/cff6dff6202e73a08c541f612940ee6b4f7b58cc))
- Switch to GH actions ([`455638e`](https://github.com/szaghi/FLAP/commit/455638eedfee696e09c36db8313bb6e77f37eea0))

## [v1.0.5](https://github.com/szaghi/FLAP/tree/v1.0.5) (2019-04-12)
[Full Changelog](https://github.com/szaghi/FLAP/compare/v1.0.4...v1.0.5)
### Miscellaneous
- Update vecfor ([`e553c4c`](https://github.com/szaghi/FLAP/commit/e553c4cda0c83075cdfa65e85ceda58c0f633953))
- Merge tag 'v1.0.4' into develop

Distance-connectivity made efficient by AABB tree

Stable release, not fully backward compatible. ([`03dfdb2`](https://github.com/szaghi/FLAP/commit/03dfdb222cf93c72ae6281611caaa3726e55c76d))
- Update vecfor ([`065a857`](https://github.com/szaghi/FLAP/commit/065a8576aeac493f228995e7a0157dbb41d09ffc))
- Fix issue[#2](https://github.com/szaghi/FLAP/issues/2)

Translate method of surface object did not recompute metrix, now fixed. ([`91e2a80`](https://github.com/szaghi/FLAP/commit/91e2a80ccaad541ddcce22f5b4b1e25513cc8f08))
- Merge branch 'release/1.0.5' ([`3b48721`](https://github.com/szaghi/FLAP/commit/3b48721053154dbefd8954affc98ba68905e2e8a))

## [v1.0.4](https://github.com/szaghi/FLAP/tree/v1.0.4) (2018-06-04)
[Full Changelog](https://github.com/szaghi/FLAP/compare/v1.0.3...v1.0.4)
### Miscellaneous
- Merge tag 'v1.0.3' into develop

Add merge files method, stable release ([`b120619`](https://github.com/szaghi/FLAP/commit/b12061947565c9e7c4e2dcf92aee89e6564dbea6))
- Temporary stop this feature ([`c7c94f0`](https://github.com/szaghi/FLAP/commit/c7c94f0ffc10a3f2038d8d0992dcc106b5f8b84f))
- Greatly improved nearby search speed by means of AABB tree ([`e2c5c45`](https://github.com/szaghi/FLAP/commit/e2c5c45e9a52dd20da45cdf15bd4744951bab2bc))
- Split STL file handler and STL surface: a more sane approach! Break backward compatibility ([`433ed02`](https://github.com/szaghi/FLAP/commit/433ed023b6e18a8fcdee8fb5c61cb55b95da718a))
- Trim out pointer procedures for not introducing too much F2008 features yet not widely supported ([`88def9d`](https://github.com/szaghi/FLAP/commit/88def9db3e3e094420afbbf59d08fe1033ccda58))
- (re)add AABB STL output ([`92ef922`](https://github.com/szaghi/FLAP/commit/92ef922f98e08319d5103ed4a3f50c7280915d84))
- Add compute_distance method that returns also facet index of closest facet ([`7d066ff`](https://github.com/szaghi/FLAP/commit/7d066ff2c32252f2502ff04d475a499d863967fe))
- Make distance/connectivity computation efficient

Make distance/connectivity computation efficient: in order to load and
analize very big input STL files the distance and connectivity
computations have strongly refactored in order to fully exploit AABB
tree representation.

The backward compatibility has been brocken: file handler is now
a distinguished object with respect the surface object.

A new method `compute_distance` for computing the distance has been added:
it returns also the facet closest to the given point, not only the distance. ([`0c92731`](https://github.com/szaghi/FLAP/commit/0c927312edf42ce2d4c843e5a63c4f8c3096bec2))
- Merge branch 'release/1.0.4' ([`cce3c21`](https://github.com/szaghi/FLAP/commit/cce3c217f9b23ac798e877dca0021afdb481b28d))

## [v1.0.3](https://github.com/szaghi/FLAP/tree/v1.0.3) (2018-05-17)
[Full Changelog](https://github.com/szaghi/FLAP/compare/v1.0.2...v1.0.3)
### Miscellaneous
- Merge tag 'v1.0.2' into develop

Add clip method to STL file handler

Add clip method to STL file handler: triangulated surface can be clipped
by means of a bounding box (axis-aligned). Optionally, it is also
possible to retain the remainder part of the surface into a separate
file handler. ([`b5573d3`](https://github.com/szaghi/FLAP/commit/b5573d38a7778db74b2116420826549a24cba4e1))
- Add method to merge STL files

Add method to merge STL files ([`f098382`](https://github.com/szaghi/FLAP/commit/f098382f01190da0cafcba3e0d0c36beb9d29dad))
- Add clip and merge to fossilizer ([`f63d63c`](https://github.com/szaghi/FLAP/commit/f63d63c15f435d1dedb5ab1441b66cfef639d901))
- Update README ([`c9c83c9`](https://github.com/szaghi/FLAP/commit/c9c83c9e5454873702fa4a2f54c02fe866a398d8))
- Merge branch 'release/1.0.3' ([`c6f6b58`](https://github.com/szaghi/FLAP/commit/c6f6b583194d2119ef951bb6a74c6e96da11b720))

## [v1.0.2](https://github.com/szaghi/FLAP/tree/v1.0.2) (2018-05-17)
[Full Changelog](https://github.com/szaghi/FLAP/compare/v1.0.1...v1.0.2)
### Miscellaneous
- Merge tag 'v1.0.1' into develop

Add standalone manipulator and automatic nearby facets reconnection

Stable release, partially backward compatible.

Add standalone manipulator (src/app/fossilizer.f90) that does quite the
same job of admesh (filling holes is only feature missing) and add the
automatic facets reconnection of nearby edges being disconnected. ([`e783ede`](https://github.com/szaghi/FLAP/commit/e783edee8816bb12a6142199d0d5d9d7027e0ee0))
- Add clip method to STL file handler

Add clip method to STL file handler: triangulated surface can be clipped
by means of a bounding box (axis-aligned). Optionally, it is also
possible to retain the remainder part of the surface into a separate
file handler. ([`02ff421`](https://github.com/szaghi/FLAP/commit/02ff421588cc9f737cac8598a0321be1711c4024))
- Merge branch 'release/1.0.2' ([`dae3ea6`](https://github.com/szaghi/FLAP/commit/dae3ea6a5c704832ac404f8415018382e8241067))

## [v1.0.1](https://github.com/szaghi/FLAP/tree/v1.0.1) (2018-05-16)
[Full Changelog](https://github.com/szaghi/FLAP/compare/v1.0.0...v1.0.1)
### Bug fixes
- Fix readme broken links ([`1d117b2`](https://github.com/szaghi/FLAP/commit/1d117b2cf15f31d0ef1b130c5e5197e0d2a78e21))

### Miscellaneous
- Merge tag 'v1.0.0' into develop

First stable release

First stable release, features:

* [X] User-friendly methods for IO STL files:
    * [x] input:
        * [x] automatic guessing of file format (ASCII or BINARY);
        * [x] load STL file effortless;
    * [x] output:
        * [x] save STL file effortless;
* [x] powerful surface analysis and manipulation:
    * [x] build facets connectivity;
    * [x] sanitize normals:
        * [x] reverse normals:
        * [x] make normals consistent:
    * [x] compute volume;
    * [x] rotate facets;
    * [x] translate facets;
    * [x] mirror facets;
    * [x] resize (scale) facets;
    * [x] compute minimal distance:
        * [x] square distance;
        * [x] square root distance;
        * [x] signed distance:
            * [x] by means of solid angle computation;
            * [x] by means of rays intersection count;
        * [x] AABB (Axis-Aligned Bounding Box) tree acceleration with user defined refinement levels;
    * [x] point-in-polyhedra test:
        * [x] by means of solid angle computation;
        * [x] by means of rays intersection count; ([`b116c2e`](https://github.com/szaghi/FLAP/commit/b116c2e5b73b59681e563b22315fc6e2bf685cf0))
- Prepare for API documentation ([`0cf10d5`](https://github.com/szaghi/FLAP/commit/0cf10d544964b71c316dbe29279dbd02c2feaa1d))
- Prepare for API documentation ([`b187948`](https://github.com/szaghi/FLAP/commit/b187948eca48e86378d0eb43689b109abfe8ee75))
- Correct README hyperlink to GH pages ([`5b2ae6c`](https://github.com/szaghi/FLAP/commit/5b2ae6cf36709cef0261d4eea37d1fafd6cafad7))
- Update building scripts ([`e2e143b`](https://github.com/szaghi/FLAP/commit/e2e143b61f5e7836fd06b0b9044f56d33c5ef561))
- Merge branch 'master' into develop ([`583d052`](https://github.com/szaghi/FLAP/commit/583d05294d991f753af790d3a9ac1e76abd7ef0a))
- Add fossilizer, a (partial) emulation of admesh ([`a2f90cf`](https://github.com/szaghi/FLAP/commit/a2f90cf992e1e5b15bb5762acb6c5d9d5b741dbc))
- Merge branch 'release/1.0.1' ([`aa0b4b0`](https://github.com/szaghi/FLAP/commit/aa0b4b04c2b647fa8a338054fbe61da5e1edea84))

## [v1.0.0](https://github.com/szaghi/FLAP/tree/v1.0.0) (2018-05-09)
[Full Changelog](https://github.com/szaghi/FLAP/compare/v0.9.5...v1.0.0)
### Miscellaneous
- Merge tag 'v0.9.5' into develop

Add repair methods

Stable release, fully backward compatible.

FOSSIL is now able to:

+ reconstruct STL connectivity
+ compute STL volume
+ make normals consistent ([`bbf6807`](https://github.com/szaghi/FLAP/commit/bbf68078dbf0c039b1cd522817c7b4e8a4a0391a))
- Add transformation methods

Add transformation methods to:

+ rotate
+ translate
+ mirror
+ resize ([`3d7fba8`](https://github.com/szaghi/FLAP/commit/3d7fba8ce235fb4c1f23ed8ca699c80b30c63394))
- Update main_page.md ([`7c48c37`](https://github.com/szaghi/FLAP/commit/7c48c373600abd3342cfc9d76ce755f42563cd1f))
- First stable release

First stable release, features:

* [X] User-friendly methods for IO STL files:
    * [x] input:
        * [x] automatic guessing of file format (ASCII or BINARY);
        * [x] load STL file effortless;
    * [x] output:
        * [x] save STL file effortless;
* [x] powerful surface analysis and manipulation:
    * [x] build facets connectivity;
    * [x] sanitize normals:
        * [x] reverse normals:
        * [x] make normals consistent:
    * [x] compute volume;
    * [x] rotate facets;
    * [x] translate facets;
    * [x] mirror facets;
    * [x] resize (scale) facets;
    * [x] compute minimal distance:
        * [x] square distance;
        * [x] square root distance;
        * [x] signed distance:
            * [x] by means of solid angle computation;
            * [x] by means of rays intersection count;
        * [x] AABB (Axis-Aligned Bounding Box) tree acceleration with user defined refinement levels;
    * [x] point-in-polyhedra test:
        * [x] by means of solid angle computation;
        * [x] by means of rays intersection count; ([`c39daf0`](https://github.com/szaghi/FLAP/commit/c39daf0af904c57535efaeaf0ef36fbc8b9f58ac))
- Merge branch 'release/1.0.0' ([`e2a110a`](https://github.com/szaghi/FLAP/commit/e2a110a109f3a7b14db00ec4cdaf1efcd63b9fd9))

## [v0.9.5](https://github.com/szaghi/FLAP/tree/v0.9.5) (2018-05-07)
[Full Changelog](https://github.com/szaghi/FLAP/compare/v0.9.0...v0.9.5)
### Miscellaneous
- Merge tag 'v0.9.0' into develop

First stable release

+ IO work;
+ basic sanitize work;
+ signed distance computation (with also AABB acceleration) work; ([`3586fcb`](https://github.com/szaghi/FLAP/commit/3586fcb3ac33ef501e2f8ec412f4f58a6cd85970))
- Add utils module ([`75c4fe2`](https://github.com/szaghi/FLAP/commit/75c4fe2625af3b9a6c3788226412be69fa09062f))
- Merge branch 'master' into develop ([`ec9be4c`](https://github.com/szaghi/FLAP/commit/ec9be4c6c2a044eac98457aabcbcd64a3146dcee))
- Add connectivity reconstruction: sanitize normal must be completed ([`435c0be`](https://github.com/szaghi/FLAP/commit/435c0be9be97568131b049eea5643a65d0d0a6c5))
- Sanitize normals work

Sanitize normals work ([`a569de1`](https://github.com/szaghi/FLAP/commit/a569de1864c3af0bdb9a4f27ff2556f5629bb008))
- Merge branch 'feature/try-to-add-stl-cleanup' into develop ([`0ab0d1c`](https://github.com/szaghi/FLAP/commit/0ab0d1c41d19c5acd6c09b993c21ce91e1c57e1b))
- Merge branch 'release/0.9.5' ([`a5b2fa2`](https://github.com/szaghi/FLAP/commit/a5b2fa208bc8be44c8b7e4b427e65a32a50c6034))

## [v0.9.0](https://github.com/szaghi/FLAP/tree/v0.9.0) (2018-05-02)
[Full Changelog](https://github.com/szaghi/FLAP/compare/v0.0.1...v0.9.0)
### Bug fixes
- Fix travis submodule issue ([`11d44e1`](https://github.com/szaghi/FLAP/commit/11d44e1bcaeda59c2c98536cf1d0c6c46dd4a48d))

### Miscellaneous
- Merge tag 'v0.0.1' into develop

First beta release. ([`7c8268f`](https://github.com/szaghi/FLAP/commit/7c8268f6843de6770348168da85869deffd9671e))
- Update travis deploy ([`75c296e`](https://github.com/szaghi/FLAP/commit/75c296e90a0e075fd08b95b328b45dce920439d4))
- Merge branch 'master' into develop ([`4e08ea7`](https://github.com/szaghi/FLAP/commit/4e08ea75a7ae8f5c7f04d0e48ab6e2384c7fb00d))
- Add normal check methods

Add methods to check and sanitize normals consistency with vertices
data. ([`3e020bd`](https://github.com/szaghi/FLAP/commit/3e020bd852e4018ea4e91ed9a8795adee912876b))
- Add distance computation: projection point region analysis must be introduced ([`190db77`](https://github.com/szaghi/FLAP/commit/190db77635735b7f96e6a23b970a12bb65c6f366))
- Signed distance computation works, must be made more efficient ([`0b3f199`](https://github.com/szaghi/FLAP/commit/0b3f199724e827cf028ca8cf595c2fe40b1f679e))
- Merge branch 'feature/add-compute-distance-method' into develop ([`190a54a`](https://github.com/szaghi/FLAP/commit/190a54a396924d3750fbd1f3e404a00fb5c1e13c))
- Add aabb_object, methods must be still implemented ([`37eefd6`](https://github.com/szaghi/FLAP/commit/37eefd6d1d553507be232dd5581fa7dd3ddd9242))
- Improve aabb object ([`390be2c`](https://github.com/szaghi/FLAP/commit/390be2c4c6e363cf1e0b0a8ec35e62e7d865c44f))
- Step over on AABB tree construction ([`8b2ce56`](https://github.com/szaghi/FLAP/commit/8b2ce567b16a99abba84be90eb752710f378bf40))
- Big step toward AABB acceleration ([`831d281`](https://github.com/szaghi/FLAP/commit/831d281d6f2739d0a651bf93cc889038710d75b8))
- Add working AABB-based distance computation

Add working AABB-based distance computation: the distance now can be
computed by brute-force search over all facets or by AABB-refinement
search.

The sing computation is still available only by brute force search. ([`bd3cc01`](https://github.com/szaghi/FLAP/commit/bd3cc01c3ac2044f344f6c719918d448ef1eeaa0))
- Update readme ([`a041129`](https://github.com/szaghi/FLAP/commit/a04112906fe82541926db5a643c66d65587dd11b))
- Add aabb extents update in initialization ([`fd5e82a`](https://github.com/szaghi/FLAP/commit/fd5e82a3fc9c500fadb822c683fb50a782b63e70))
- Remove not useful aabb_cloud search ([`043864e`](https://github.com/szaghi/FLAP/commit/043864ee60c92264ae323513893b6b088674e500))
- Add full AABB acceleration

Add full AABB acceleration: both distance and sign computations can now
be accelerated by means of AABB tree refinement. ([`7c84e42`](https://github.com/szaghi/FLAP/commit/7c84e423242ca4ab974c29ed328d482e148a4c7b))
- Merge branch 'feature/add-aabb-ray-intersection-test' into develop ([`38b1cb6`](https://github.com/szaghi/FLAP/commit/38b1cb6390a33d6a656dfd057c4ab92876f66776))
- Update fobos coverage ([`a2fb5af`](https://github.com/szaghi/FLAP/commit/a2fb5af37da485d3c17a1d0f2cd242c08df30985))
- Merge branch 'release/0.9.0' ([`2165b69`](https://github.com/szaghi/FLAP/commit/2165b69ffc7c5b772066e04066b4ae3855ad262e))

## [v0.0.1](https://github.com/szaghi/FLAP/tree/v0.0.1) (2018-04-10)
### Miscellaneous
- Init

First commit of FOSSIL library, a pure Fortran library to parse STL
files. ([`56d5266`](https://github.com/szaghi/FLAP/commit/56d5266c89001e63f6c1c00ba3fe5579f67f928b))
- Add facet object and much more

Add facet object and much more:
+ add facet object
+ add file stl object
+ add file stl IO with test in both ascii and binary forms ([`870ce82`](https://github.com/szaghi/FLAP/commit/870ce8228ce4f582a42b96edcd5178a9da85f97f))
- FOSSIL first beta

FOSSIL first beta: Fortran parser of STL files. ([`24cb162`](https://github.com/szaghi/FLAP/commit/24cb1621f30bdaf8f93b449c8f8c7edba07d9ff0))
- Merge branch 'release/0.0.1' ([`3cc41f4`](https://github.com/szaghi/FLAP/commit/3cc41f45a8d6c6fe81b611af432cfdff92da8923))


