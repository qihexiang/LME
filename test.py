from typing import Any
from layers.StaticLayer import StaticLayer
from layers.UtilLayers import AutoBondLayer
from workflow.Workflow import Workflow, default_runners
from fs.osfs import OSFS
import yaml
from libs.RadiiTable import default_radius_table

workflow_data = yaml.load(OSFS("./").readtext("test.yml"), yaml.Loader)

class ClampRisky:
    def __init__(self, options, meta):
        self.radius_table = {
            key: default_radius_table[key] for key in default_radius_table
        }
    
    def __call__(self, item, tags):
        re_connected = StaticLayer(contains=(item, [AutoBondLayer()]))
        if len(re_connected.bond_ids) != len(item.bond_ids):
            print(f"risky: tags: {tags}")
            return item, "warning"
        return item, ""     

workflow = Workflow(workflow_data, runners=default_runners | {"risky": (ClampRisky, False)})

workflow.run()