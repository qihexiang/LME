import yaml
from fs.osfs import OSFS
from os.path import isabs, join
from layers.EditableLayer import EditableLayer
from pydash import py_

class Workflow:
    def __init__(self, config, runners) -> None:
        version = config["version"]
        rootDirectory = config["rootDirectory"]
        subsitutes = config["substitutes"]
        subsitutes = [item if isabs(item) else join(self.rootDirectory, item) for item in subsitutes]
        self.metas = {
            "version": version, "rootDirectory": rootDirectory,
            "subsitutes": subsitutes
        }
        self.template = config["template"]
        self.template = self.template if isabs(self.template) else join(self.rootDirectory, self.template)
        self.template = EditableLayer.from_mol2(OSFS("/").readtext(self.template)).to_static_layer()
        self.generated = [(self.template, {})]
        self.jobs = config["jobs"]
        self.current_step = 0
        self.runners = runners
    
    def run_step(self):
        if(self.current_step >= len(self.jobs)):
            raise None
        current_task = self.jobs[self.current_step]
        self.current_step += 1
        (runner_builder, need_flat) = self.runners[current_task["use"]]
        runner = runner_builder(current_task["with"], self.metas)
        def full_runner(working_item):
            model, names = working_item
            if need_flat:
                products = runner(model, names)
                products = [(product, names | {current_task["name"]: tag}) for product, tag in products]
            else:
                model, tag = runner(model, names)
                return (model, names | {current_task["name"]: tag})
        processed = py_.map(self.generated, full_runner)
        if need_flat:
            processed = py_.flatten(processed)
        self.generated = processed
        return None

    def run(self):
        while(True):
            try:
                self.run_step()
            except None:
                break
            except Exception as err:
                print(err)
                exit(1)
