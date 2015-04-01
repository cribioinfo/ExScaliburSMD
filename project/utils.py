'''Utility functions for the ExScalibur-SMD pipeline'''
from collections import OrderedDict
import os
import yaml

def check_dir(path):
    '''Checks for/Creates a path.

    -- path - The path to check/create
    '''
    if not os.path.isdir(path):
        os.makedirs(path)

def check_file(fil, fname):
    '''Checks if a file can be read. Throws an exception.'''
    try:
        fp = open(fil)
    except IOError as e:
        raise IOError('Unable to access {0}: {1}\n'.format(fname, fil))
 
def represent_odict(dump, tag, mapping, flow_style=None):
    """For writing out OrderedDicts for easier user-end configuration.

    This is not my original code and was provided by Michael Elsdorfer
    https://gist.github.com/miracle2k/3184458#file-odict-py
    """
    value = []
    node = yaml.MappingNode(tag, value, flow_style=flow_style)
    if dump.alias_key is not None:
        dump.represented_objects[dump.alias_key] = node
    best_style = True
    if hasattr(mapping, 'items'):
        mapping = mapping.items()
    for item_key, item_value in mapping:
        node_key = dump.represent_data(item_key)
        node_value = dump.represent_data(item_value)
        if not (isinstance(node_key, yaml.ScalarNode) and not node_key.style):
            best_style = False
        if not (isinstance(node_value, yaml.ScalarNode) and not node_value.style):
            best_style = False
        value.append((node_key, node_value))
    if flow_style is None:
        if dump.default_flow_style is not None:
            node.flow_style = dump.default_flow_style
        else:
            node.flow_style = best_style
    return node

yaml.SafeDumper.add_representer(OrderedDict,
    lambda dumper, value: represent_odict(dumper, u'tag:yaml.org,2002:map', value))
