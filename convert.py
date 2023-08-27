from fs.osfs import OSFS
from pydash import py_
from libs.RadiiTable import RadiiTable

number_to_element = lambda x: RadiiTable[RadiiTable["Z"]==x]["Element"].values[0]

file = OSFS(".").readtext("data.xyz")
lines = file.split("\r\n")
lines = py_.filter(lines, lambda line: line != "")
new_content = f"@<TRIPOS>MOLECULE\r\nname\r\n{len(lines)} {0}\r\nSMALL\r\nNO_CHARGES\r\n\r\n@<TRIPOS>ATOM\r\n"
for i, line in enumerate(lines):
    [ele, x, y, z] = py_.filter(line.split(" "), lambda item: item != "")
    ele = number_to_element(int(ele))
    new_content += f"{i+1} {ele}{i+1} {x} {y} {z} {ele}\r\n"

OSFS(".").writetext("converted.mol2", new_content)

