#include "exceptions.h"

namespace combi_ff {


xml_error::xml_error(const std::string& xml_tag, const std::string& message) :
	input_error("XML Error: in XmlElement '" + xml_tag + "' -> " + message) {}

opening_tag_error::opening_tag_error(const std::string& filename, const std::string& xml_tag, const std::string& word) :
	xml_error(xml_tag, "opening_tag_error: in file " + filename + " expected opening tag '<', but instead got " + word)  {}

closing_tag_error::closing_tag_error(const std::string& filename, const std::string& xml_tag, const std::string& word) :
	xml_error(xml_tag, "closing_tag_error: in file " + filename + " expected closing tag '>', but instead got " + word)  {}

mismatched_quotes_error::mismatched_quotes_error(const std::string& filename, const std::string& xml_tag, const std::string& word,
                                                 const char& openingQuote, const char& closingQuote) :
	xml_error(xml_tag, "mismatched_quotes_error: in file " + filename + " opening quote is " + std::string(1, openingQuote)
	          + " but closing quote is " + std::string(1, closingQuote) + " in " + word) {}

empty_tag_error::empty_tag_error(const std::string& filename, const std::string& xml_tag, const std::string& word) :
	xml_error(xml_tag, "empty_tag_error: in file " + filename + " encountered empty tag: " + word + "\n") {}

unexpected_word_error::unexpected_word_error(const std::string& filename, const std::string& xml_tag, const std::string& expected,
                                             const std::string& got) :
	xml_error(xml_tag, "unexpected_word_error: in file " + filename + " expected " + expected + " but got " + got + "\n") {}

duplicate_attributes_error::duplicate_attributes_error(const std::string& filename, const std::string& xml_tag, std::string& attributeType,
                                                       std::string& attributeValuePrev, std::string& attributeValueCur) :
	xml_error(xml_tag, "duplicate_attributes_error: in file " + filename + " encountered duplicate attribute type " + attributeType + " (once with " + attributeValuePrev +
	          " and once with " + attributeValueCur + "\n") {}

xml_version_error::xml_version_error(const std::string& filename, const std::string& xml_tag, const std::string& versionString,
                                     const char ch) :
	xml_error(xml_tag, "xml_version_error: in file " + filename +
	          " version number has to be a number, but got " + versionString + " (invalid character '" + ch + "')") {}

root_tag_mismatch_error::root_tag_mismatch_error(const std::string& filename, std::string& found, std::string& expected) :
	xml_error(found, "root_tag_mismatch_error: in file " + filename + " expected root element to be " + expected +
	          " according to !DOCTYPE info, but got " + found + "\n") {}

unexpected_closing_tag_error::unexpected_closing_tag_error(const std::string& filename, std::string& found,
                                                           const std::string expected) :
	xml_error(found, "unexpected_closing_tag_error: in file " + filename + " expected " + expected + " but got " + found + "\n") {}

wrong_root_tag_error::wrong_root_tag_error(const std::string& expected, const std::string& found) :
	xml_error(found, "wrong_root_tag_error: expected root tag " + expected + " but got " + found + "\n") {}

unexpected_end_of_file_error::unexpected_end_of_file_error(const std::string& filename, const std::string& xml_tag,
                                                           const std::string& lastwords) :
	xml_error(xml_tag, "unexpected_end_of_file_error: in file " + filename + " unexpected end of file at: '" + lastwords + "'\n") {}

incomplete_comment_error::incomplete_comment_error(const std::string& filename, const std::string& word) :
	xml_error("comment", "incomplete_comment_error: in file " + filename + " encountered unexpected '!' in " + word + " (for comments, please use '<!-- ... -->')\n") {}

unexpected_tag_error::unexpected_tag_error(const std::string& found, const std::string& expected) :
	xml_error(found, "unexpected_tag_error: encountered unexpected tag '" + found + " (expected '" + expected + "')\n") {}

missing_attribute_error::missing_attribute_error(const std::string& filename, const std::string& tag,
                                                 const std::string& attribute) :
	xml_error(tag, "missing_attribute_error: in file " + filename + " missing attribute " + attribute + " in element " + tag + "\n") {}

wrong_amount_of_children_error::wrong_amount_of_children_error(const std::string& tag, const std::string& expected,
                                                               const size_t found) :
	xml_error(tag, "wrong_amount_of_children_error: wrong number of children in element " + tag + " (expected " + expected + ", but got " + std::to_string(
	              found) + ")\n") {}

missing_value_error::missing_value_error(const std::string& xml_tag) :
	xml_error(xml_tag, "missing_value_error: missing value in element " + xml_tag + "\n") {}

wrong_amount_of_attributes_error::wrong_amount_of_attributes_error(const std::string& xml_tag, const size_t expected,
                                                                   const size_t found) :
	xml_error(xml_tag, "wrong_amount_of_attributes_error: wrong number of attributes in element " + xml_tag + " (expected "
	          + std::to_string(expected) + " but found " + std::to_string(found) + ")\n") {}

unexpected_value_error::unexpected_value_error(const std::string& xml_tag, const std::string& value) :
	xml_error(xml_tag, "unexpected_value_error: found value " + value + " for element " +
	          xml_tag + " but this element is expected to have no value\n") {}




} //namespace combi_ff