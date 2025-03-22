workspace(name = "openroad")

# TODO: move that to a use_repo_rule() in MODULE.bazel, and/or upstream
# to bazel

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

http_archive(
    name = "tk_tcl",
    build_file = "//bazel/tk_tcl:bundled.BUILD.bazel",
    sha256 = "5196dbf6638e3df8d5c87b5815c8c2b758496eb6f0e41446596c9a4e638d87ed",
    strip_prefix = "tcl8.6.10",
    urls = [
        "https://prdownloads.sourceforge.net/tcl/tcl8.6.10-src.tar.gz",
    ],
)
