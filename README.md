# Tristan-MP v2.1.1

For detailed tutorials and code description please visit our [wiki](https://ntoles.github.io/tristan-wiki/). If you are a new user testing the code on a new computer cluster, please consider contributing to [this chapter](https://ntoles.github.io/tristan-wiki/tristanv2-configure.html#cluster-specific-customization) about cluster-specific configuration to make the life easier for future generations.

### Branch policy

To prevent `v2` from growing to become the Lovecraftian monster it once was we highly encourage both users and developers to follow the guidlines on branching politics.

* For development:
    * all the new features shall be added to `dev/<feature>` branch;
    * as soon as the new features are tested they can be pushed to the main development branch: `dev/main`.
* For users:
    * user-specific branches are allowed (e.g. to test userfiles), but as soon as it works we highly encourage people to merge their new userfiles to `dev`;
    * the naming for user-specific branches shall be the following: `user/<problem>` or `user/<username>`;
    * if any of the core routines is modified, we highly encourage to use the `dev/...` branching instead of the `user/...`.

The `dev/main` branch is the most up-to-date **tested** version of the code, and it will be merged to `master` as soon as all the new features are documented.

### Latest updates

We employ [semantic versioning](https://semver.org/) for this code. Given a version number `v<MAJOR>.<MINOR>.<PATCH>`, increment the:
- `MAJOR` version when you make incompatible changes,
- `MINOR` version when you add functionality in a backwards compatible manner, and
- `PATCH` version when you make backwards compatible bug fixes.

* `v2.0.1` __Jul 2020__
  * Compton scattering and GCA pusher added
  * Slice outputs for 3d added
  * Particle momentum binning improved for downsampling
  * Minor bugs fixed
* `v2.1` __Dec 2020__
  * Pair annihilation module (+ advanced test for all QED modules coupled)
  * Merging of charged particles
  * Particle payloads
  * Vay pusher
  * Individual particle current deposition (necessary for reflecting walls and downsampling of charged particles)
  * Spatial binning for spectra
  * Automated testing framework with fulltest.py
  * Major restructuring of output modules
  * Major restructuring of QED modules
  * Minor issues fixed for GCA pusher
* `v2.1.1` __Jan 2021__
  * Patched the Compton cross section normalization
  * Minor bug fixes
