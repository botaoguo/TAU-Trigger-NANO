def LoadAsList(file_name):
    trig_vpset = []
    tag_path_names = []
    import json
    with open(file_name) as f:
        trig_desc = json.load(f)
    for trig_name, desc in trig_desc.iteritems():
        filters = [ str(','.join(path_list)) for path_list in desc['filters'] ]
        is_tag = 'is_tag' in desc and desc['is_tag'] > 0
        leg_types = [ str(leg_type) for leg_type in desc['leg_types'] ]
        path = str(trig_name)
        #print("PATH : {}".format(path))
        trig_vpset.append( (filters, is_tag, leg_types, path) )
        if is_tag:
            tag_path_names.append( str(trig_name) )
    return trig_vpset, tag_path_names

class TriggerDescriptor:
    
    def __init__(self):
        self.name = 'triggerdescriptor'
        
class TriggerLeg:
    
    def __init__(self):
        self.name = "triggerleg"
        
class TriggerResult:
    
    def __init__(self):
        self.name = "trigresult"