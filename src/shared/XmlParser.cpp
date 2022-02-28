// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "XmlParser.h"

#include <fstream>

namespace combi_ff {

/*********************
XmlElement definitions
*********************/

XmlElement::XmlElement(XmlElement* parent)
    : parent(parent),
      distance(parent->distance + 1),
      file_name(parent->file_name) {}
XmlElement::XmlElement(XmlElement* parent, const std::string& tag)
    : parent(parent),
      distance(parent->distance + 1),
      file_name(parent->file_name),
      tag(tag) {}
XmlElement::XmlElement(XmlElement* parent, const XmlElement& element)
    : parent(parent),
      distance(parent->distance + 1),
      file_name(parent->file_name),
      tag(element.tag),
      attributes(element.attributes),
      value(element.value) {}

void XmlElement::AddElement() {
  children.push_back(XmlElement_ptr(new XmlElement(this)));
}
void XmlElement::AddElement(const XmlElement& element) {
  children.push_back(XmlElement_ptr(new XmlElement(this, element)));

  // std::cout << "  " << children.back()->tag << std::endl;
  for (auto&& child : element.children) children.back()->AddElement(*child);
}
void XmlElement::AddElement(const std::string& tag) {
  children.push_back(XmlElement_ptr(new XmlElement(this, tag)));
}
void XmlElement::AddElement(const std::string& tag,
                            const std::list<XmlElement_ptr>::iterator pos) {
  children.insert(pos, XmlElement_ptr(new XmlElement(this, tag)));
}
void XmlElement::CheckTagName(const std::string& tag_) const {
  if (tag != tag_) throw unexpected_tag_error(tag, tag_);
}
void XmlElement::CheckAttributeSize(const size_t n) const {
  if (attributes.size() != n)
    throw wrong_amount_of_attributes_error(tag, n, attributes.size());
}
void XmlElement::CheckAttribute(const std::string& attribute) const {
  if (attributes.find(attribute) == attributes.end())
    throw missing_attribute_error(*file_name, tag, attribute);
}
void XmlElement::CheckNumberOfChildren_equal(const size_t n) const {
  if (children.size() != n)
    throw wrong_amount_of_children_error(tag, std::to_string(n),
                                         children.size());
}
void XmlElement::CheckNumberOfChildren_atLeast(const size_t n) const {
  if (children.size() < n)
    throw wrong_amount_of_children_error(tag, ">=" + std::to_string(n),
                                         children.size());
}
void XmlElement::CheckNoValue() const {
  if (value.size()) throw unexpected_value_error(tag, value);
}
void XmlElement::CheckValue() const {
  if (!value.size()) throw missing_value_error(tag);
}

void XmlElement::SetAttribute(const XmlElement& el,
                              const std::string& attributeName) {
  auto&& attribute = el.attributes.find(attributeName);
  if (attribute != el.attributes.end())
    attributes[attributeName] = attribute->second;
}

void XmlElement::SetFileName(const std::string& file_name_) {
  file_name = std::shared_ptr<const std::string>(new std::string(file_name_));
}

XmlElement_ptr XmlElement::GetLastChild() const { return children.back(); }

XmlElement_ptr XmlElement::GetFirstChild() const { return children.front(); }

void XmlElement::RemoveLastChild() {
  auto&& toBeErased = children.end();
  children.erase(--toBeErased);
}

void XmlElement::RemoveChildren() { children.clear(); }

void XmlElement::RemoveChild(std::list<XmlElement_ptr>::iterator child) {
  children.erase(child);
}

const size_t XmlElement::GetNumberOfChildren() const { return children.size(); }

/******************
XmlTree definitions
******************/

XmlTree::XmlTree(const std::string& roottag) : root() { root.tag = roottag; }
XmlElement& XmlTree::GetRoot() { return root; }
const XmlElement& XmlTree::GetRoot_const() const { return root; }
XmlElement* XmlTree::GetRootPointer() { return &root; }
void XmlTree::SetRoot(const std::string& root_tag,
                      const std::map<std::string, std::string>& attributes) {
  root.tag = root_tag;
  root.attributes = attributes;
}
void XmlTree::CheckRootTagName(const std::string& root_tag) const {
  if (root.tag != root_tag) throw wrong_root_tag_error(root_tag, root.tag);
}
void XmlTree::AddSubtree(const XmlTree& tree) {
  root.AddElement(tree.GetRoot_const());
}

/********************
XmlParser definitions
********************/

XmlParser::XmlParser(const std::string& file_name) : file_name(file_name) {
  tree.GetRoot().SetFileName(file_name);
}
XmlParser::XmlParser(const std::string& file_name, const std::string& root_tag_)
    : tree(root_tag_),
      file_name(file_name),
      root_tag(root_tag_),
      dtd(root_tag + ".dtd") {
  tree.GetRoot().SetFileName(file_name);
}

XmlTree& XmlParser::GetTree() { return tree; }
XmlElement& XmlParser::GetRoot() { return tree.GetRoot(); }

const XmlTree& XmlParser::GetTree_const() const { return tree; }
const std::string& XmlParser::GetVersion() const { return version; }
const std::string& XmlParser::GetEncoding() const { return encoding; }
const std::string& XmlParser::GetRootTag() const { return root_tag; }
const std::string& XmlParser::GetDtd() const { return dtd; }

/**********************
XmlParserIn definitions
**********************/

XmlParserIn::XmlParserIn(const std::string& file_name, const xml_read_mode mode)
    : XmlParser(file_name), xml_file(file_name), read_mode(mode) {
  if (!xml_file.is_open())
    throw std::runtime_error("XML input file " + file_name + " is not open\n");

  if (read_mode == read_all)
    CreateTree();

  else {
    ReadXMLInfo();
    ReadDocInfo();
    CheckRootTag();
    cur_node = tree.GetRootPointer();
  }
}
XmlParserIn::~XmlParserIn() { xml_file.close(); }

void XmlParserIn::CreateTree() {
  ReadXMLInfo();
  ReadDocInfo();
  BuildTree();
}

void XmlParserIn::ReadDocInfo() {
  auto start_pos = xml_file.tellg();
  std::string tag_content = ReadNextTag();

  if (tag_content.size() < 2 ||
      tag_content[1] != '!') {  // we know that size of tag_content is at least
                                // 3 (checking in readNextTag)
    std::cout << "?Warning: didn't find any DOCTYPE info\n";
    xml_file.seekg(start_pos);
    return;
  }

  std::string word = ReadNextWord(tag_content);

  if (word != "<!DOCTYPE")
    throw unexpected_word_error(file_name, tag_content, "<!DOCTYPE", word);

  root_tag = ReadNextWord(tag_content);

  if (root_tag.back() == '>') {
    root_tag = root_tag.substr(0, root_tag.size() - 1);
    std::cout << "?Warning: no DTD file given\n";
  }

  else {
    word = ReadNextWord(tag_content);

    if (word != "SYSTEM")
      throw unexpected_word_error(file_name, tag_content, "SYSTEM", word);

    dtd = ReadNextWord(tag_content);

    if (dtd.back() == '>') {
      dtd = dtd.substr(0, dtd.size() - 1);
      tag_content = ">";
    }

    char opening_quote = dtd.front();
    char closing_quote = dtd.back();

    if ((opening_quote != '\'' && opening_quote != '\"') ||
        opening_quote != closing_quote)
      throw mismatched_quotes_error(file_name, "tag_content", dtd,
                                    opening_quote, closing_quote);

    dtd = dtd.substr(1, dtd.size() - 2);

    if (tag_content != ">")
      throw closing_tag_error(file_name, tag_content, tag_content);
  }
}

void XmlParserIn::BuildTree() {
  CheckRootTag();
  cur_node = tree.GetRootPointer();

  if (!xml_file.eof()) {
    std::string tag_content = ReadNextTag();

    while (tag_content != "</" + tree.GetRoot().tag + ">") {
      // std::cout << "cur_node: " << cur_node->tag << std::endl;
      // std::cout << "distance: " << cur_node->distance << std::endl;
      if (tag_content.size() > 1 && tag_content[1] == '/') {
        if (tag_content.substr(2, cur_node->tag.size()) != cur_node->tag)
          throw unexpected_closing_tag_error(file_name, tag_content,
                                             "</" + cur_node->tag + ">");

        cur_node = cur_node->parent;
        tag_content = ReadNextTag();
      }

      else if (tag_content == "/>") {
        cur_node = cur_node->parent;
        tag_content = ReadNextTag();
      }

      else {
        cur_node->AddElement();
        cur_node = &*cur_node->children.back();
        std::string elementName = ReadElementName(tag_content);
        cur_node->tag = elementName;

        if (tag_content != ">") {
          std::map<std::string, std::string> attributes =
              ReadAttributes(elementName, tag_content);
          cur_node->attributes = attributes;
        }

        if (tag_content != "/>") {
          std::string value = ReadElementValue();
          tag_content = ReadNextTag();
          cur_node->value = value;
        }
      }
    }
  }
}

std::string XmlParserIn::ReadElementValue() {
  std::string value("");

  while (xml_file.peek() != '<') {
    if (xml_file.eof())
      throw unexpected_end_of_file_error(file_name, cur_node->tag, value);

    // if(xml_file.peek() != ' ' && xml_file.peek() != '\t')
    value += (char)xml_file.get();
    // else xml_file.Get();
  }

  auto pos = xml_file.tellg();
  xml_file.get();

  // ignore comments
  if (xml_file.peek() == '!') {
    // std::cout << "ignore comment\n";
    std::string comment("");
    // read in !--
    comment += (char)xml_file.get();
    comment += (char)xml_file.get();
    comment += (char)xml_file.get();

    if (comment == "!--") {
      char ch(' ');

      while (ch != '>' && comment.substr(comment.size() - 3, 3) != "-->") {
        xml_file >> ch;
        comment += ch;
      }
    }

    else
      throw incomplete_comment_error(file_name, value + comment);

    value += ReadElementValue();
  }

  else if (xml_file.peek() != '/') {
    value = "";
    xml_file.seekg(pos);
  }

  else
    xml_file.seekg(pos);

  if (value.size() && value.back() == '\n')
    value = value.substr(0, value.size() - 1);

  return value;
}

void XmlParserIn::CheckRootTag() {
  std::string tag_content = ReadNextTag();
  std::string element_name = ReadElementName(tag_content);
  std::map<std::string, std::string> attributes;

  if (tag_content != ">")
    attributes = ReadAttributes(element_name, tag_content);

  if (element_name != root_tag && root_tag != "")
    throw root_tag_mismatch_error(file_name, element_name, root_tag);

  tree.SetRoot(element_name, attributes);

  if (tag_content == "/>") std::cout << "only root tag found\n";
}

std::map<std::string, std::string> XmlParserIn::ReadAttributes(
    const std::string& tag, std::string& tag_content) {
  std::map<std::string, std::string> attributes;

  while (tag_content != ">" && tag_content != "/>") {
    std::string attribute_type = ReadAttributeType(tag_content);
    std::string attribute_value = ReadAttributeValue(tag_content);

    if (attributes.find(attribute_type) != attributes.end())
      throw duplicate_attributes_error(file_name, tag, attribute_type,
                                       attribute_value,
                                       attributes[attribute_type]);

    else
      attributes.emplace(attribute_type, attribute_value);
  }

  return attributes;
}
std::string XmlParserIn::ReadNextTag() {
  char ch;
  std::string tag_content;
  xml_file >> ch;

  if (ch != '<') {
    xml_file >> tag_content;
    throw opening_tag_error(file_name, std::string(1, ch) + tag_content,
                            std::string(1, ch) + tag_content);
  }

  else {
    tag_content = ch;

    while (!xml_file.eof() && xml_file.peek() != '>') {
      ch = (char)xml_file.get();
      tag_content += ch;

      if (ch == ' ' || ch == '\t') {
        while (!xml_file.eof() &&
               (xml_file.peek() == ' ' || xml_file.peek() == '\t'))
          xml_file.get();
      }
    }

    if (xml_file.eof())
      throw closing_tag_error(file_name, tag_content, tag_content);

    else {  // case that xml_file.peek() == '>'
      xml_file >> ch;
      tag_content += ch;
    }
  }

  if (tag_content.size() < 3)
    throw empty_tag_error(file_name, tag_content, tag_content);

  // ignore comments
  if (tag_content.substr(0, 4) == "<!--") {
    while (!xml_file.eof() &&
           tag_content.substr(tag_content.size() - 3, 3) != "-->") {
      while (xml_file.peek() != '>') {
        xml_file >> ch;
        tag_content += ch;
      }

      xml_file >> ch;
      tag_content += ch;
    }

    if (xml_file.eof())
      throw unexpected_end_of_file_error(file_name, tag_content, tag_content);

    return ReadNextTag();
  }

  return tag_content;
}
std::string XmlParserIn::ReadNextWord(std::string& tag_content) {
  std::string word("");
  size_t i = 0;

  while (tag_content[i] == ' ' ||
         tag_content[i] == '\t')  // not really necessary, because readNextTag
                                  // removes duplicate whitespaces
    i++;

  for (; i < tag_content.size(); i++) {
    if (tag_content[i] != ' ' && tag_content[i] != '\t')
      word += tag_content[i];

    else {
      tag_content = tag_content.substr(i + 1, tag_content.size());
      return word;
    }
  }

  if (word.back() != '>')
    throw closing_tag_error(file_name, tag_content, tag_content);

  return word;
}
std::string XmlParserIn::ReadElementName(std::string& tag_content) {
  if (tag_content.front() != '<')
    throw opening_tag_error(file_name, tag_content, tag_content);

  std::string element_name = ReadNextWord(tag_content);

  if (element_name ==
      "<")  // case that there's a whitespace between '<' and element name
    element_name = ReadNextWord(tag_content);

  else
    element_name = element_name.substr(1, element_name.size());

  if (element_name.back() == '>') {
    if (element_name.size() > 1 && *(element_name.end() - 2) == '/') {
      tag_content = "/>";
      element_name = element_name.substr(0, element_name.size() - 2);
    }

    else {
      tag_content = ">";
      element_name = element_name.substr(0, element_name.size() - 1);
    }
  }

  return element_name;
}
std::string XmlParserIn::ReadAttribute(const std::string& attribute_keyword,
                                       std::string& tag_content,
                                       const TagType type) {
  std::string word = ReadNextWord(tag_content);

  if (word.size() < attribute_keyword.size() ||
      word.substr(0, attribute_keyword.size()) != attribute_keyword)
    throw unexpected_word_error(file_name, tag_content,
                                attribute_keyword + "= (attribute)", word);

  else if (word.size() == attribute_keyword.size()) {
    word = ReadNextWord(tag_content);

    if (word.front() != '=')
      throw unexpected_word_error(file_name, tag_content, "=", word);

    if (word.size() == 1)
      word = ReadNextWord(tag_content);

    else
      word = word.substr(1, word.size());
  }

  else {
    if (word.back() == '=')
      word = ReadNextWord(tag_content);

    else {
      if (word[attribute_keyword.size()] != '=')
        throw unexpected_word_error(file_name, tag_content, "=", word);

      word = word.substr(attribute_keyword.size() + 1, word.size());
    }
  }

  char opening_quote = word.front();
  char closing_quote;

  if (word.back() == '>')
    closing_quote = (type == normal ? *(word.end() - 2) : *(word.end() - 3));

  else
    closing_quote = word.back();

  if ((opening_quote != '\"' && opening_quote != '\'') ||
      closing_quote != opening_quote)
    throw mismatched_quotes_error(file_name, tag_content, word, opening_quote,
                                  closing_quote);

  std::string attribute;

  if (word.back() == '>') {
    if (type == normal) {
      attribute = word.substr(1, word.size() - 3);
      tag_content = ">";
    }

    else {
      attribute = word.substr(1, word.size() - 4);
      tag_content = word.substr(word.size() - 2, word.size());
    }
  }

  else
    attribute = word.substr(1, word.size() - 2);

  return attribute;
}
std::string XmlParserIn::ReadAttributeType(std::string& tag_content) {
  std::string attribute_type = ReadNextWord(tag_content);
  std::size_t pos = attribute_type.find('=');

  if (pos == std::string::npos) {
    std::string word = ReadNextWord(tag_content);

    if (word.front() != '=')
      throw unexpected_word_error(file_name, tag_content, "=",
                                  word + " " + tag_content);

    if (word.size() > 1)
      tag_content = word.substr(1, word.size()) + " " + tag_content;
  }

  else {
    tag_content = attribute_type.substr(pos + 1, attribute_type.size()) + " " +
                  tag_content;
    attribute_type = attribute_type.substr(0, pos);
  }

  return attribute_type;
}
std::string XmlParserIn::ReadAttributeValue(std::string& tag_content) {
  std::string attribute_value = ReadNextWord(tag_content);

  if (attribute_value.back() == '>') {
    if (*(attribute_value.end() - 2) == '/') {
      tag_content = "/>";
      attribute_value = attribute_value.substr(0, attribute_value.size() - 2);
    }

    else {
      tag_content = ">";
      attribute_value = attribute_value.substr(0, attribute_value.size() - 1);
    }
  }

  char opening_quote = attribute_value.front();
  char closing_quote = attribute_value.back();

  if ((opening_quote != '\'' && opening_quote != '\"') ||
      opening_quote != closing_quote)
    throw mismatched_quotes_error(file_name, tag_content, attribute_value,
                                  opening_quote, closing_quote);

  attribute_value = attribute_value.substr(1, attribute_value.size() - 2);
  return attribute_value;
}
void XmlParserIn::ReadXMLInfo() {
  std::string tag_content = ReadNextTag();

  if (tag_content[1] != '?') {  // we know that size of tag_content is at least
                                // 3 (checking in readNextTag)
    std::cout
        << "?Warning: didn't find any info on version (using default version "
        << version << ") or encoding (using default encoding " << encoding
        << ")\n";
    xml_file.seekg(0);
    return;
  }

  std::string word = ReadNextWord(tag_content);

  if (word != "<?xml")
    throw unexpected_word_error(file_name, tag_content, "<?xml", word);

  std::string version_string =
      ReadAttribute("version", tag_content, question_mark);
  bool hasDecimalPoint(false);

  for (auto&& ch : version_string) {
    if (!std::isdigit(ch)) {
      if (ch == '.' && !hasDecimalPoint)
        hasDecimalPoint = true;

      else
        throw xml_version_error(file_name, tag_content, version_string, ch);
    }
  }

  if (tag_content == "?>") {
    std::cout << "?Warning: didn't find encoding info, using default "
              << encoding << "\n";
    return;
  }

  encoding = ReadAttribute("encoding", tag_content, question_mark);
  word = ReadNextWord(tag_content);

  if (word != "?>")
    throw unexpected_word_error(file_name, tag_content, "?> (closing tag)",
                                word);
}

bool XmlParserIn::ReadUntilElement(const std::string& read_until_tag_name) {
  if (cur_node->children.size()) RemoveChildrenOfCurrentNode();

  if (!xml_file.eof()) {
    std::string tag_content = ReadNextTag();
    if (tag_content == "</" + cur_node->tag + ">") {
      return false;
    }

    while (tag_content != "</" + read_until_tag_name + ">") {
      // std::cout << "cur_node: " << cur_node->tag << std::endl;
      // std::cout << "distance: " << cur_node->distance << std::endl;
      if (xml_file.eof())
        throw unexpected_end_of_file_error(file_name, tag_content, tag_content);

      if (tag_content.size() > 1 && tag_content[1] == '/') {
        if (tag_content.substr(2, cur_node->tag.size()) != cur_node->tag)
          throw unexpected_closing_tag_error(file_name, tag_content,
                                             "</" + cur_node->tag + ">");

        cur_node = cur_node->parent;
        tag_content = ReadNextTag();
        //				std::cout << "    a " << tag_content <<
        // std::endl;

      }

      else if (tag_content == "/>") {
        cur_node = cur_node->parent;
        tag_content = ReadNextTag();
        // std::cout << "    b " << tag_content << std::endl;
      }

      else {
        cur_node->AddElement();
        cur_node = &*cur_node->children.back();
        std::string element_name = ReadElementName(tag_content);
        cur_node->tag = element_name;

        if (tag_content != ">") {
          std::map<std::string, std::string> attributes =
              ReadAttributes(element_name, tag_content);
          cur_node->attributes = attributes;
        }

        if (tag_content != "/>") {
          std::string value = ReadElementValue();
          tag_content = ReadNextTag();
          // std::cout << "    c " << tag_content << std::endl;
          cur_node->value = value;
        }
      }
    }

    cur_node = cur_node->parent;
    return true;
  }

  return false;
}

void XmlParserIn::RemoveLastChildOfCurrentNode() {
  cur_node->RemoveLastChild();
}

void XmlParserIn::RemoveChildrenOfCurrentNode() { cur_node->RemoveChildren(); }

bool XmlParserIn::eof() const { return xml_file.eof(); }

/***********************
XmlParserOut definitions
***********************/

XmlParserOut::XmlParserOut(const std::string& file_name,
                           const std::string& root_tag)
    : XmlParser(file_name, root_tag), xml_file(file_name) {
  if (!xml_file.is_open())
    throw std::runtime_error("XML output file " + file_name + " is not open\n");
}
XmlParserOut::~XmlParserOut() { xml_file.close(); }

void XmlParserOut::WriteHead() {
  xml_file << "<?xml version='" << version << "' encoding='" << encoding
           << "'?>\n"
           << "<!DOCTYPE " << root_tag << " SYSTEM \"" << dtd << "\">\n"
           << "<" << tree.GetRoot().tag;

  if (tree.GetRoot().attributes.size()) {
    for (auto&& at : tree.GetRoot().attributes)
      xml_file << " " << at.first << "=\"" << at.second << "\"";
  }

  xml_file << ">\n";
}

void XmlParserOut::WriteTail() {
  xml_file << "</" << tree.GetRoot().tag << ">";
}

void XmlParserOut::WriteAndRemoveLastChild() {
  xml_file << *(tree.GetRoot().GetLastChild());
  tree.GetRoot().RemoveLastChild();
}

/*****************
ofstream operators
*****************/
std::ostream& operator<<(std::ostream& stream, const XmlElement& el) {
  std::string whitespaces(el.distance * xml_indent_size, ' ');
  stream << whitespaces << "<" << el.tag;

  if (el.attributes.size()) {
    for (auto&& at : el.attributes)
      stream << " " << at.first << "=\"" << at.second << "\"";
  }

  // if(!el.children.size() && !el.value.size())
  // stream << "/>\n";

  // else {
  if (!el.children.size())
    stream << ">" << el.value << "</" << el.tag << ">\n";

  else {
    stream << ">\n";

    if (el.value.size()) stream << whitespaces << "  " << el.value << "\n";

    for (auto&& ch : el.children) stream << *ch;

    stream << whitespaces << "</" << el.tag << ">\n";
    //}
  }

  return stream;
}

std::ostream& operator<<(std::ostream& stream, const XmlTree& tree) {
  stream << tree.GetRoot_const() << std::endl;
  return stream;
}

std::ostream& operator<<(std::ostream& stream, const XmlParser& parser) {
  stream << "<?xml version='" << parser.GetVersion() << "' encoding='"
         << parser.GetEncoding() << "'?>\n";
  stream << "<!DOCTYPE " << parser.GetRootTag() << " SYSTEM \""
         << parser.GetDtd() << "\">\n";
  stream << parser.GetTree_const() << std::endl;
  return stream;
}

}  // namespace combi_ff