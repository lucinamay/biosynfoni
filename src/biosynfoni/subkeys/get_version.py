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

from biosynfoni.subkeys.biosmartfonis import substructureSmarts
from biosynfoni.subkeys.versionfonis import fpVersions
from biosynfoni.subkeys.default_version import defaultVersion


def get_values(
    value_name: str,
    version: str = defaultVersion,
    subs_smarts: dict = substructureSmarts,
    fp_versions: dict[str, list[str]] = fpVersions,
) -> list[str]:
    """gives list of values of choice of all substructures in fingerprint version
    input:   fp_version_name (str) -- name of the version
             subs_smarts (dict) {substructure ids (e.g. 'fp1') :  {keys,values}}
             fp_versions (dict) {version name (e.g. fps_full_2): (list) substructure names (e.g. 'fp1')}
    output:  list of the requested values for each substructure in the version
    """
    values = []
    sub_ids = fp_versions[version]
    for sub_id in sub_ids:
        sub_info = subs_smarts[sub_id]
        if value_name not in sub_info:
            raise ValueError(
                f"Value {value_name} not found in substructure {sub_id} of version {version}"
            )
        else:
            values.append(sub_info[value_name])

    assert len(values) == len(sub_ids)
    return values


get_smarts = partial(get_values, value_name="smarts")
get_names = partial(get_values, value_name="name")
get_pathway = partial(get_values, value_name="pathway")
