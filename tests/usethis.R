# install.packages("devtools")
# install.packages("usethis")


View(iris)

# install.packages("devtools")
# install.packages("usethis")
devtools::clean_dll()
devtools::document()
devtools::install()
devtools::load_all()

usethis::use_git()
usethis::create_package(path = "~/web/ACMGuru")
devtools::load_all()
devtools::build()
devtools::install()
devtools::load_all()
devtools::check()

# Build the manual
package_path <- "~/web/ACMGuru"
devtools::build_manual(pkg = package_path, path = file.path(package_path, "./inst/doc"))


# updates ----
#  [semantic versioning](https://semver.org/) principles (major.minor.patch).

# Update the `NAMESPACE` and the documentation files in the `man/` directory:
# devtools::document()

# Update the Version Number
# Edit the `DESCRIPTION` file
# - **Patch** version for backward-compatible bug fixes.
# - **Minor** version for new features that are backward-compatible.
# - **Major** version for changes that break backward compatibility.

# **Rebuild the Package**:
# This command creates a tarball of your package that can be installed.
# devtools::build()

# **Check the Package**:
# devtools::check()

# **Install the Package Locally**:
# devtools::install()

# **Release the Package**:
#  [CRAN submission process](https://cran.r-project.org/submit.html). For a GitHub release, push your changes to GitHub and consider using GitHub releases to tag the version and provide release notes.

# **Tag the Release**:
# git tag -a v1.0.0 -m "Release version 1.0.0"
# git push --tags

# Remember, the process of updating the version number is manual to give you control over versioning semantics.

