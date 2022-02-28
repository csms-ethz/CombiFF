// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef XmlParser_H_
#define XmlParser_H_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>

#include "exceptions.h"

namespace combi_ff {

struct XmlElement;
class XmlTree;
class XmlParser;
class XmlParserIn;
class XmlParserOut;

typedef enum { normal, question_mark } TagType;
struct XmlElement;
typedef std::shared_ptr<XmlElement> XmlElement_ptr;

std::ostream& operator<<(std::ostream& stream, const XmlElement& el);
std::ostream& operator<<(std::ostream& stream, const XmlTree& tree);
std::ostream& operator<<(std::ostream& stream, const XmlParser& parser);

static const std::size_t xml_indent_size(2);
static const std::string xml_indent(xml_indent_size, ' ');

struct XmlElement : std::enable_shared_from_this<XmlElement> {
  XmlElement() = default;
  explicit XmlElement(XmlElement* parent);
  XmlElement(XmlElement* parent, const std::string& tag);
  XmlElement(XmlElement* parent, const XmlElement& element);

  void AddElement();
  void AddElement(const XmlElement& element);
  void AddElement(const std::string& tag);
  void AddElement(const std::string& tag,
                  const std::list<XmlElement_ptr>::iterator pos);
  void CheckTagName(const std::string& tag_) const;
  void CheckAttributeSize(const size_t n) const;
  void CheckAttribute(const std::string& attribute) const;
  void CheckNumberOfChildren_equal(const size_t n) const;
  void CheckNumberOfChildren_atLeast(const size_t n) const;
  void CheckNoValue() const;
  void CheckValue() const;
  void SetAttribute(const XmlElement& el, const std::string& attribute_name);
  void SetFileName(const std::string& file_name);
  XmlElement_ptr GetLastChild() const;
  XmlElement_ptr GetFirstChild() const;
  void RemoveLastChild();
  void RemoveChildren();
  void RemoveChild(std::list<XmlElement_ptr>::iterator child);
  const size_t GetNumberOfChildren() const;

  XmlElement* parent{NULL};
  std::list<XmlElement_ptr> children{std::list<XmlElement_ptr>(0)};
  size_t distance{0};
  std::shared_ptr<const std::string> file_name{NULL};
  std::string tag{""};
  std::map<std::string, std::string> attributes{};
  std::string value{""};
};

class XmlTree {
 public:
  XmlTree() = default;
  explicit XmlTree(const std::string& root_tag);
  XmlElement& GetRoot();
  const XmlElement& GetRoot_const() const;
  XmlElement* GetRootPointer();
  void SetRoot(const std::string& root_tag_,
               const std::map<std::string, std::string>& attributes);
  void CheckRootTagName(const std::string& root_tag) const;
  void AddSubtree(const XmlTree& tree);

 private:
  XmlElement root{};
};

class XmlParser {
 public:
  XmlParser() = default;
  explicit XmlParser(const std::string& file_name);
  XmlParser(const std::string& file_name, const std::string& root_tag);

  XmlTree& GetTree();
  XmlElement& GetRoot();

  const XmlTree& GetTree_const() const;
  const std::string& GetVersion() const;
  const std::string& GetEncoding() const;
  const std::string& GetRootTag() const;
  const std::string& GetDtd() const;

 protected:
  XmlTree tree;
  std::string version{"1.0"};     // default version is 1.0
  std::string encoding{"UTF-8"};  // default encoding is UTF-8
  std::string file_name{""};
  std::string root_tag{""};  // string for doctype (= root element)
  std::string dtd{""};       // string for dtd file
};

class XmlParserIn : public XmlParser {
 public:
  typedef enum { read_all, read_element_by_element } xml_read_mode;
  XmlParserIn(const std::string& file_name, const xml_read_mode mode);
  ~XmlParserIn();

  void CreateTree();
  void ReadXMLInfo();
  void ReadDocInfo();
  void BuildTree();
  void CheckRootTag();
  std::map<std::string, std::string> ReadAttributes(const std::string& tag,
                                                    std::string& tag_content);
  std::string ReadNextTag();
  std::string ReadNextWord(std::string& tag_content);
  std::string ReadAttribute(const std::string& attribute_keyword,
                            std::string& tag_content, const TagType type);
  std::string ReadElementName(std::string& tag_content);
  std::string ReadAttributeType(std::string& tag_content);
  std::string ReadAttributeValue(std::string& tag_content);
  std::string ReadElementValue();
  bool ReadUntilElement(const std::string& read_until_tag_name);
  void RemoveLastChildOfCurrentNode();
  void RemoveChildrenOfCurrentNode();
  const XmlElement* GetCurNode() { return cur_node; }
  bool eof() const;

 private:
  std::ifstream xml_file;
  xml_read_mode read_mode{read_all};
  XmlElement* cur_node{NULL};
};

class XmlParserOut : public XmlParser {
 public:
  XmlParserOut(const std::string& file_name, const std::string& root_tag);
  ~XmlParserOut();

  void WriteHead();
  void WriteTail();
  void WriteAndRemoveLastChild();

 private:
  std::ofstream xml_file;
};

}  // namespace combi_ff

#endif