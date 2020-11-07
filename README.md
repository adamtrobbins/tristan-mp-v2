# Tristan-MP v2

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
