from workflow.Workflow import Workflow
from fs.osfs import OSFS
import yaml

workflow_data = yaml.load(OSFS("./").readtext("test.yml"), yaml.Loader)

workflow = Workflow(workflow_data)

workflow.run()