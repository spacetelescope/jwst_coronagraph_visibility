# Release Procedure

The `jwst_coronagraph_visibility` team employs the following procedure for creating a new release:

1. Update appropriate version numbers in locations that store the version number
2. Update the release notes
3. Open, review, and merge pull request with the release procedure changes
4. Create a new tag/release on GitHub
5. Upload new version of software to PyPI
6. Upload new version of software to astroconda

Detailed instructions for performing a release are given below:

### 1. Update the version number in in locations that store the version number

Update the version number in all locations of hard-coded version numbers. For example, the VERSION variable in `setup.py` 
should be updated to the new version number of the release, using the [x.y.z convention](#version-conv).

### 3. Update the release notes

In [CHANGES](CHANGES.md), write a concise but detailed description of all of the notable changes that have
occurred since the last release. One way to acquire this information is to scroll through the commit history of
the project, and look for commits in which a pull request was merged.


### 4. Open, review, and merge pull requests with the release procedure changes

Once you've committed the changes from (1), (2), and (3) in your branch, push your branch to GitHub/GitLab using
the upstream remote, open a pull request and assign reviewers. Either you or the reviewer should eventually merge these pull
requests.

### 5. Create a new tag/release on GitHub/GitLab

Once the pull request into the production branch from (4) has been merged, click on the releases button on the
main page of the repository, then hit the "Draft a new release button". The "Tag version" should be the version
number of the release, the "Target" should be the production branch, the "Release title" should (also) be the
version number of the release, and the "Description" should match that of the changelog entry in (3). Once all
of that information is added, hit the big green "Publish" release button.

### 6. Upload new version of software to PyPI

 - Follow the `instructions from PyPA <https://packaging.python.org/distributing/#uploading-your-project-to-pypi>`_ to set up ``twine`` and login credentials.
 - Check out the just-tagged commit: ``git checkout 0.x.y``
 - Create the sdist: ``python setup.py sdist``
 - Create a universal wheel: ``python setup.py bdist_wheel --universal``
 - Check distribution is valid: ``twine check dist/*``
 - Test upload with ``twine upload --repository-url https://test.pypi.org/legacy/ dist/*``
 - Upload to PyPI: ``twine upload dist/*``


### 7. Upload new version of software to astroconda

 - Fork `astroconda/astroconda-contrib <https://github.com/astroconda/astroconda-contrib>`_ to your GitHub, if you haven't already.
 - Clone your fork locally (or pull upstream changes, if you already have a local clone).
 - Check out a feature branch like this: ``git checkout -b jwst_coronagraph_visibility-0.x.y``.
 - Update the version number in ``astroconda-contrib/jwst_coronagraph_visibility/meta.yaml`` (along with any version requirements on upstream dependencies).
 - Run `conda build -c http://ssb.stsci.edu/astroconda --skip-existing --python=3.7 jwst_coronagraph_visibility-0.x.y`
  and verify it passes without error
 - Commit your changes.
 - Push to GitHub.
 - Open a pull request against `astroconda/astroconda-contrib <https://github.com/astroconda/astroconda-contrib>`_.


