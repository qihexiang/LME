import yaml
from fs.osfs import OSFS
from posixpath import isabs, join, relpath
from datetime import datetime
from layers.EditableLayer import EditableLayer
from pydash import py_
from workflow.runners import default_runners

class Workflow:
    def __init__(self, config, runners = default_runners) -> None:
        version = config["version"]
        rootDirectory = join(*config["rootDirectory"])
        subsitutes = py_.map(config["substitutes"], lambda paths: join(*paths))
        subsitutes = [relpath(item) if isabs(item) else item for item in subsitutes]
        self.metas = {
            "version": version, "rootDirectory": rootDirectory,
            "substitutes": subsitutes, "output_prefix": config["output_prefix"]
        }
        template_path = join(*config["template"])
        template_path = relpath(template_path, rootDirectory) if isabs(template_path) else template_path
        self.rootDirectory = OSFS(rootDirectory)
        self.template = EditableLayer.from_mol2(self.rootDirectory.readtext(template_path)).to_static_layer()
        self.generated = [(self.template, {})]
        self.jobs = config["jobs"]
        self.current_step = 0
        self.runners = runners
    
    def run_step(self):
        if(self.current_step >= len(self.jobs)):
            return False
        current_task = self.jobs[self.current_step]
        self.current_step += 1
        print(f"Enter job {self.current_step}: {current_task['name']}")
        start_at = datetime.now()
        (runner_builder, need_flat) = self.runners[current_task["use"]]
        runner = runner_builder(current_task["with"], self.metas)
        def full_runner(working_item):
            model, names = working_item
            if need_flat:
                products = runner(model, names)
                products = [(product, names | {current_task["name"]: tag}) for product, tag in products]
                return products
            else:
                model, tag = runner(model, names)
                return (model, names | {current_task["name"]: tag})
        processed = py_.map(self.generated, full_runner)
        if need_flat:
            processed = py_.flatten(processed)
        self.generated = processed
        print(f"Exit job {self.current_step}: {current_task['name']}, uses {datetime.now() - start_at}")
        return True

    def run(self):
        start_at = datetime.now()
        print(f"Task start at {start_at}")
        while(True):
            if not self.run_step():
                break
        print(f"All jobs finished. Used {datetime.now() - start_at}")
