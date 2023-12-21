"""
||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: biosynfoni arrangement       ||
created: 2023-09                    ||
author: Lucina-May Nollen           ||
institute: WUR Bioinformatics       ||
____________________________________

||||||||||||  ()()()  |||||||||||||||

Biosynfoni Definition. Can be converted into an sdf file per 'version'

"""
from functools import partial

from biosynfoni.biosmartfonis import Substructures as SUBSTRUCTURES
from biosynfoni.versionfonis import fpVersions as FP_VERSIONS
from biosynfoni.versionfonis import defaultVersion

DEFAULT_BIOSYNFONI_VERSION = defaultVersion  # used as default in concerto-fp


def get_values(
    fp_version_name: str,
    value_name: str,
    subs_smarts: dict = SUBSTRUCTURES,
    fp_versions: dict[str, list[str]] = FP_VERSIONS,
) -> list[str]:
    """gives list of values of choice of all substructures in fingerprint version
    input:   fp_version_name (str) -- name of the version
             subs_smarts (dict) {substructure ids (e.g. 'fp1') :  {keys,values}}
             fp_versions (dict) {version name (e.g. fps_full_2): (list) substructure names (e.g. 'fp1')}
    output:  list of the requested values for each substructure in the version
    """
    chosen_sub_names = fp_versions[fp_version_name]
    return [[x, subs_smarts[x][value_name]] for x in chosen_sub_names]


get_smarts = partial(get_values, value_name="smarts")
get_names = partial(get_values, value_name="name")

# def get_smarts(
#     fp_version_name: str,
#     subs_smarts: dict = SUBSTRUCTURES,
#     fp_versions: dict[str, list[str]] = FP_VERSIONS,
# ) -> list[list[str, str]]:
#     """gives list of smarts of substructures of choice
#     input:   fp_version_name (str) -- name of the version
#              subs_smarts (dict)
#                          (key) substructure names (e.g. 'fp1')
#                          (val) (dict)
#                          'smarts': substructure RDK molfiles (f/SMARTS)
#              fp_versions (dict)
#                          (key) version name (e.g. fps_full_2)
#                          (val) (list) substructure names (e.g. 'fp1')
#     output:  list of the requested values for each substructure in the version
#     """
#     return

# def get_names(
#     fp_version_name: str,
#     subs_smarts: dict = SUBSTRUCTURES,
#     fp_versions: dict[str, list[str]] = FP_VERSIONS,
# ) -> list[str]:
#     """gives list of names of substructures of choice
#     input:   fp_version_name (str) -- name of the version
#              subs_smarts (dict)
#                          (key) substructure names (e.g. 'fp1')
#                          (val) (dict)
#                          'smarts': substructure RDK molfiles (f/SMARTS)
#              fp_versions (dict)
#                          (key) version name (e.g. fps_full_2)
#                          (val) (list) substructure names (e.g. 'fp1')
#     """
#     chosen_sub_names = fp_versions[fp_version_name]
#     return [[x, subs_smarts[x]["name"]] for x in chosen_sub_names]


def main():
    # get_smarts(DEFAULT_BIOSYNFONI_VERSION, SUBSTRUCTURES, FP_VERSIONS)
    # get_subsset(DEFAULT_BIOSYNFONI_VERSION, SUBSTRUCTURES, FP_VERSIONS)

    # get_subsset('regular_1007')
    return


if __name__ == "__main__":
    main()
