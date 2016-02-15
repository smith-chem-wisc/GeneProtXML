# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "anthony"
__date__ = "$Oct 28, 2015 6:50:01 PM$"

def set_prefixes(root, prefix_map):
    def fixup_element_prefixes(element, uri_map, memo):
        def fixup(name):
            try:
                return memo[name]
            except KeyError:
                if name[0] != "{":
                    return
                uri, tag = name[1:].split("}")
                if uri in uri_map:
                    new_name = uri_map[uri] + ":" + tag
                    memo[name] = new_name
                    return new_name
            #fix element name
            name = fixup(element.tag)
            if name:
                element.tag = name
            #fix attribute names
            for key, value in element.items():
                name = fixup(key)
                if name:
                    element.set(name, value)
                    del element.attrib[key]
    
    #build uri map and add to root element
    uri_map = {}
    for prefix, uri in prefix_map.items():
        uri_map[uri] = prefix
        root.set("xmlns:" + prefix, uri)
        
    #fixup all elements in teh tree
    memo = {}
    for element in root.getiterator():
        fixup_element_prefixes(element, uri_map, memo)