# -*- Python -*-

"""Bazel rule to run flex toolchain
"""

# Adapter rule around the @rules_flex toolchain.
def genlex(name, prefix, src, out):
    """Generate C/C++ language source from lex file using Flex
    """
    # how can I pass PREFIX_OPTION to the cmd ?
    PREFIX_OPTION='-P ' + prefix if prefix else ''
    if prefix:
        print("TODO(hzeller): add flex prefix support {}".format(PREFIX_OPTION))
    native.genrule(
        name = name,
        srcs = [src],
        outs = [out],
        cmd = select({
            "//bazel:use_local_flex_bison_enabled": "flex --outfile=$@  $<",
            "@platforms//os:windows": "win_flex.exe --outfile=$@ $<",
            "//conditions:default": "M4=$(M4) $(FLEX) --outfile=$@  $<",
        }),
        toolchains = select({
            "//bazel:use_local_flex_bison_enabled": [],
            "@platforms//os:windows": [],
            "//conditions:default": [
                "@rules_flex//flex:current_flex_toolchain",
                "@rules_m4//m4:current_m4_toolchain",
            ],
        }),
    )
