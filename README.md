# Tristan-MP v2

For a detailed tutorials and descriptions please visit our [wiki](https://ntoles.github.io/tristan-wiki/).

### Branch politics

To prevent `v2` from growing to become the Lovecraftian monster it once was we highly encourage both users and developers to follow the guidlines on branching politics. 

* For development:
    * all the new features shall be added to `dev/<feature>` branch;
    * as soon as the new features are tested they can be pushed to the main development branch: `dev`.
* For users:
    * user-specific branches are allowed (e.g. to test userfiles), but as soon as it works we highly encourage people to merge their new userfiles to `dev`;
    * the naming for user-specific branches shall be the following: `user/<problem>` or `user/<username>`;
    * if any of the core routines is modified, we highly encourage to user `dev/...` branch instead of `user/...`.

The `dev` branch is the most up-to-date **tested** version of the code, and it will be merged to `master` as soon as all the new features are documented.
