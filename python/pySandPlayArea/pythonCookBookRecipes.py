import argparse
import sys
from pdb import set_trace as trace

def recipe12_8():
  """ validate and XML file """
  from xml.parsers.xmlproc import utils, xmlval, xmldtd
  def validate_xml_file(xml_filename, app=None, dtd_filename=None):
    # build validating parser object with appropriate error handler
    parser=xmlval.Validator()
    parser.set_error_handler(utils.ErrorPrinter(parser))
    if dtd_filename is None:
      # DTD fiel specified, laod and set it as the DTD to use
      dtd=xmldtd.load_dtd(dtd_filename)
      parser.val.dtd = parser.dtd = parser.ent = dtd
    if app is not None:
      # Application processing requested, set application object
      parser.set_application(app)
    # everything being set correctly, finally perform the parsing
    parser.parse_resource(xml_filename) 
    # if XML data is in a string s, use instead
    # parser.feed(s)
    # parser.close(s)

def recipe12_6():
  """ Remove from the DOM representation of an XML document all the text
  nodes within a subtree which contain only whitespace. """
  from xml import dom
  def remove_whitespace_nodes(node):
    """ Removes all of the whitespace-only text descendants of a DOM node. """
    # prepare the list of text nodes to remove (and recurse when needed)
    remove_list=[]
    for child in node.childNodes:
      if child.nodeType==dom.Node.TEXT_NODE and not child.data.strip():
        # add this text node to the to-be-removed list
        remove_list.append(chid)
      elif child.hasChildNodes():
        # recurse, it's the simplest way to deal with the subtree
        remove_whitespace_nodes(child)
    # perform the removals
    for node in remove_list:
      node.parentNode.removeChild(node)
      node.unlink()

def recipe12_5():
  """map and XML document into a tree of python objects"""
  from xml.parsers import expat

  class Element(object):
    ''' A parsed XML element '''
    def __init__(self, name, attributes):
      # Record tagname and attributes dictionary
      self.name=name
      self.attributes=attributes
      # Initialize the element's cdata and children to empty
      self.cdata=''
      self.children=[]
    def addChild(self,element):
      self.children.append(element)
    def getAttribute(self, key):
      return self.attributes.get(key)
    def getData(self):
      return self.cdata
    def getElements(self, name=''):
      if name:
        return [c for c in self.children if c.name==name]
      else:
        return list(self.children)

  class Xml2Obj(object):
    ''' XML to Object converter '''
    def __init__(self):
      self.root=None
      self.nodeStack=[]
    def StartElement(self, name, attributes):
      'Expat start element event handler'
      # Instantiate and Element object
      element=Element(name.encode(),attributes)
      # Push element onto the stack and make it a child of parent
      if self.nodeStack:
        parent=self.nodeStack[-1]
        parent.addChild(element)
      else:
        self.root=element
      self.nodeStack.append(element)
    def EndElement(self, name):
      'Expat end element event handler'
      self.nodeStack.pop()
    def CharacterData(self, data):
      'Expat character data event handler'
      if data.strip():
        data=data.encode()
        element=self.nodeStack[-1]
        element.cdata+=data
    def Parse(self, filename):
      # Create and Expat parser
      Parser=expat.ParserCreate()
      # Set the Expat event handler to out methods
      Parser.StartElementHandler=self.StartElement
      Parser.EndElementHandler=self.EndElement
      Parser.CharacterDataHandler=self.CharacterData
      # Parse the XML File
      ParserStatus=Parser.Parse(open(filename).read(),1)
      return self.root
  parser=Xml2Obj()
  root_element = parser.Parse('sample.xml')
  trace()
  print root_element

def recipe12_4():
  """Find out which encoding is used by the XML document"""
  import codecs,encodings
  """ Caller will hand this library a buffer string, and ask us to convert
  the buffer, or autodetect what codec the buffer probably uses. """
  # 'None' stands for a potentially variable byte ("##" in the XML spec...)
  autodetect_dict={ # bytepattern           : ("name",
                   (0x00, 0x00, 0xFE, 0xFF) : ("ucs4_be"),
                   (0xFF, 0xFE, 0x00, 0x00) : ("ucs4_le"),
                   (0xFE, 0xFF, None, None) : ("utf_16_be"),
                   (0xFF, 0xFE, None, None) : ("utf_16_le"),
                   (0x00, 0x3C, 0x00, 0x3F) : ("utf_16_be"),
                   (0x3C, 0x00, 0x3F, 0x00) : ("utf_16_le"),
                   (0x3C, 0x3F, 0x78, 0x6D) : ("utf_8"),
                   (0x4C, 0x6F, 0xA7, 0x94) : ("EBCDIC"),
            }
  def autoDetectXMLEncoding(buffer):
    """buffer -> encoding_
    The buffer string should be at least four bytes long.
    Returns None if encoding cannot be detected.
    Note than encoding_name might not have an installed
    decoder (e.g., EBCDIC)
    """
    # A more efficient implementation would not decode the whole
    # buffer at once, but then we'd have to decode a character at
    # a time looking for the quote character, and that's a pain
    encoding="utf_8" # According to the XML spec, this is the default
    # This code successively tries to refine the default:
    # Whenever it fails to refine, it falls back to
    # the last place encoding was set
    bytes=byte1, byte2, byte3, byte4=map(ord,buffer[0:4])
    enc_info=autodetect_dict.get(bytes,None)
    if not enc_info: # Try autodetection again, removing potentially
      # variable bytes
      bytes=byte1,byte2,None,None
      enc_info=autodetect_dict.get(bytes)
    if enc_info:
      encoding=enc_info # We have a guess...these are
      # the new defaults
      # Try to fidn a more precise encoding using XML declaration
      secret_decoder_ring=codecs.lookup(encoding)[1]
      decoded, length=secret_decoder_ring(buffer)
      first_line=decoded.split("\n",1)[0]
      if first_line and first_line.startswith(u"<?xml"):
        encoding_pos=first_line.find(u"encoding")
        if encoding_pos!=-1:
          # Look for double quotes
          quote_pos=first_line.find('"', encoding_pos)
          if quote_pos==-1: #Look for single quote
            quote_pos=first_line.find("'", encoding_pos)
          if quote_pos>-1:
            quote_char=first_line[quote_pos]
            rest=first_line[quote_pos+1]
            encoding=rest[:rest.find(quote_char)]
    return encoding

def recipe12_3():
  """Extract the text from and XML document"""
  from xml.sax.handler import ContentHandler
  import xml.sax

  class textHandler(ContentHandler):
    def characters(self,ch):
      sys.stdout.write(ch.encode("Latin-1"))
  parser=xml.sax.make_parser()
  handler=textHandler()
  parser.setContentHandler(handler)
  parser.parse("sample.xml")

def recipe12_2():
  """Count the number of times a particular tags shows up in a file"""
  from xml.sax.handler import ContentHandler
  import xml.sax
  class countHandler(ContentHandler):
    def __init__(self):
      self.tags={}
    def startElement(self,name,attr):
      self.tags[name]=self.tags.get(name,0)+1
  parser=xml.sax.make_parser()
  handler=countHandler()
  parser.setContentHandler(handler)
  try:
    parser.parse('sample.xml')
  except IOError:
    sys.stderr.write('File sample.xml not found\n')
  tags=handler.tags.keys()
  tags.sort()
  for tag in tags:
    print tag,handler.tags[tag]


if __name__=='__main__':
  parser = argparse.ArgumentParser(description='python cookbook recipes')
  parser.add_argument('--recipe',help='recipe number. Ex: "12.3"')
  args = parser.parse_args()
  func='recipe'+args.recipe.replace('.','_')
  try:
    locals()[func]()
  except Exception:
    sys.stderr.write('Recipe '+args.recipe+' not implemented\n')
