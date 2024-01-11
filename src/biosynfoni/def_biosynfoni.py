"""
||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: arrangefoni                  ||
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
from biosynfoni.default_version import defaultVersion


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


def main():
    return


if __name__ == "__main__":
    main()
