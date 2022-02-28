// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <stdexcept>

namespace combi_ff {

class input_error : public std::logic_error {
  using std::logic_error::logic_error;
};

class input_warning : public std::logic_error {
  using std::logic_error::logic_error;
};

class help_exception : public std::exception {};

class xml_error : public input_error {
 public:
  xml_error(const std::string& xml_tag, const std::string& message);
};

class opening_tag_error : public xml_error {
 public:
  opening_tag_error(const std::string& filename, const std::string& xml_tag,
                    const std::string& word);
};

class closing_tag_error : public xml_error {
 public:
  closing_tag_error(const std::string& filename, const std::string& xml_tag,
                    const std::string& word);
};

class mismatched_quotes_error : public xml_error {
 public:
  mismatched_quotes_error(const std::string& filename,
                          const std::string& xml_tag, const std::string& word,
                          const char& openingQuote, const char& closingQuote);
};

class empty_tag_error : public xml_error {
 public:
  empty_tag_error(const std::string& filename, const std::string& xml_tag,
                  const std::string& word);
};

class unexpected_word_error : public xml_error {
 public:
  unexpected_word_error(const std::string& filename, const std::string& xml_tag,
                        const std::string& expected, const std::string& got);
};

class duplicate_attributes_error : public xml_error {
 public:
  duplicate_attributes_error(const std::string& filename,
                             const std::string& xml_tag,
                             std::string& attributeType,
                             std::string& attributeValuePrev,
                             std::string& attributeValueCur);
};

class xml_version_error : public xml_error {
 public:
  xml_version_error(const std::string& filename, const std::string& xml_tag,
                    const std::string& versionString, const char ch);
};

class root_tag_mismatch_error : public xml_error {
 public:
  root_tag_mismatch_error(const std::string& filename, std::string& found,
                          std::string& expected);
};

class unexpected_closing_tag_error : public xml_error {
 public:
  unexpected_closing_tag_error(const std::string& filename, std::string& found,
                               const std::string& expected);
};

class wrong_root_tag_error : public xml_error {
 public:
  wrong_root_tag_error(const std::string& expected, const std::string& found);
};

class unexpected_end_of_file_error : public xml_error {
 public:
  unexpected_end_of_file_error(const std::string& filename,
                               const std::string& xml_tag,
                               const std::string& lastwords);
};

class incomplete_comment_error : public xml_error {
 public:
  incomplete_comment_error(const std::string& filename,
                           const std::string& word);
};

class unexpected_tag_error : public xml_error {
 public:
  unexpected_tag_error(const std::string& found, const std::string& expected);
};

class missing_attribute_error : public xml_error {
 public:
  missing_attribute_error(const std::string& filename, const std::string& tag,
                          const std::string& attribute);
};

class wrong_amount_of_children_error : public xml_error {
 public:
  wrong_amount_of_children_error(const std::string& tag,
                                 const std::string& expected,
                                 const size_t found);
};

class missing_value_error : public xml_error {
 public:
  explicit missing_value_error(const std::string& xml_tag);
};

class wrong_amount_of_attributes_error : public xml_error {
 public:
  wrong_amount_of_attributes_error(const std::string& xml_tag,
                                   const size_t expected, const size_t found);
};

class unexpected_value_error : public xml_error {
 public:
  unexpected_value_error(const std::string& xml_tag, const std::string& value);
};

}  // namespace combi_ff

#endif