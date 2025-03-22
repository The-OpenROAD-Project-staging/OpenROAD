load("@bazel_skylib//rules:common_settings.bzl", "bool_flag")

exports_files([
  "tcl_encode.bzl",
  "flex.bzl",
  "bison.bzl",
  "tcl_wrap_cc.bzl",
])

cc_library(
    name = "flex",
    deps = [
        "@rules_flex//flex:current_flex_toolchain",  # For FlexLexer.h
    ],
    visibility = ["//visibility:public"],
)

bool_flag(
    name = "use_local_flex_bison",
    build_setting_default = False,
)

config_setting(
    name = "use_local_flex_bison_enabled",
    flag_values = {":use_local_flex_bison": "true"},
)
