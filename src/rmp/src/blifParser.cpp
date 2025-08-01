// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "rmp/blifParser.h"

#include <boost/bind/bind.hpp>
#include <boost/config/warning_disable.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/container.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/define_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_alternative.hpp>
#include <iostream>
#include <string>
#include <vector>

namespace blif_parser {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

using namespace boost::placeholders;

using boost::spirit::ascii::space_type;
using boost::spirit::ascii::string;
using boost::spirit::qi::lit;
using qi::lexeme;

using ascii::char_;
using ascii::space;

void setNewInput(std::string input, rmp::BlifParser* parser)
{
  if (input != "\\") {
    parser->addInput(input);
  }
}

void setNewOutput(std::string output, rmp::BlifParser* parser)
{
  if (output != "\\") {
    parser->addOutput(output);
  }
}

void setNewClock(std::string clock, rmp::BlifParser* parser)
{
  if (clock != "\\") {
    parser->addClock(clock);
  }
}

void setNewInstanceType(std::string type, rmp::BlifParser* parser)
{
  parser->addNewInstanceType(type);
}

void setNewGate(std::string gate, rmp::BlifParser* parser)
{
  parser->addNewGate(gate);
}

void setGateNets(std::string net, rmp::BlifParser* parser)
{
  parser->addConnection(net);
}

void endParser(std::string end, rmp::BlifParser* parser)
{
  parser->endParser();
}

bool parse(std::string::iterator first,
           std::string::iterator last,
           rmp::BlifParser* parser)
{
  qi::rule<std::string::iterator, std::string(), ascii::space_type> _string;
  _string %= lexeme[+(char_ - (' ' | qi::eol))];

  qi::rule<std::string::iterator, space_type> rule
      = (lit(".model") >> _string >> lit(".inputs")
         >> +(!(&lit(".outputs"))
              >> _string[boost::bind(&setNewInput, _1, parser)])
         >> lit(".outputs")
         >> +(!(&lit(".gate")) >> !(&lit(".clock"))
              >> _string[boost::bind(&setNewOutput, _1, parser)])
         >> -(lit(".clock")
              >> +(!(&lit(".gate"))
                   >> _string[boost::bind(&setNewClock, _1, parser)]))
         >> +((lit(".gate")[boost::bind(&setNewInstanceType, "gate", parser)]
               | lit(".mlatch")[boost::bind(
                   &setNewInstanceType, "mlatch", parser)])
              >> _string[boost::bind(&setNewGate, _1, parser)]
              >> +(!(&lit(".gate")) >> !(&lit(".mlatch")) >> !(&lit(".end"))
                   >> _string[boost::bind(&setGateNets, _1, parser)]))
         >> string(".end")[boost::bind(&endParser, _1, parser)]);

  bool valid = qi::phrase_parse(first, last, rule, space);
  return valid;
}

}  // namespace blif_parser

namespace rmp {

BlifParser::BlifParser()
{
  combCount_ = 0;
  flopCount_ = 0;
  currentInstanceType_ = GateType::None;
  currentGate_ = "";
}
void BlifParser::addInput(const std::string& input)
{
  inputs_.push_back(input);
}
void BlifParser::addOutput(const std::string& output)
{
  outputs_.push_back(output);
}
void BlifParser::addClock(const std::string& clock)
{
  clocks_.push_back(clock);
}
void BlifParser::addNewInstanceType(const std::string& type)
{
  if (currentInstanceType_ != GateType::None) {
    gates_.push_back(
        Gate(currentInstanceType_, currentGate_, currentConnections_));
  }
  currentInstanceType_ = GateType::Mlatch;
  if (type != "mlatch") {
    currentInstanceType_ = (type == "gate") ? GateType::Gate : GateType::None;
  }
  if (currentInstanceType_ == GateType::Mlatch) {
    flopCount_++;
  } else if (currentInstanceType_ == GateType::Gate) {
    combCount_++;
  }
  currentConnections_.clear();
}
void BlifParser::addNewGate(const std::string& cell_name)
{
  currentGate_ = cell_name;
}
void BlifParser::addConnection(const std::string& connection)
{
  currentConnections_.push_back(connection);
}
void BlifParser::endParser()
{
  if (currentInstanceType_ != GateType::None) {
    gates_.push_back(
        Gate(currentInstanceType_, currentGate_, currentConnections_));
  }
}

const std::vector<std::string>& BlifParser::getInputs() const
{
  return inputs_;
}
const std::vector<std::string>& BlifParser::getOutputs() const
{
  return outputs_;
}
const std::vector<std::string>& BlifParser::getClocks() const
{
  return clocks_;
}
const std::vector<Gate>& BlifParser::getGates() const
{
  return gates_;
}
int BlifParser::getCombGateCount() const
{
  return combCount_;
}
int BlifParser::getFlopCount() const
{
  return flopCount_;
}

// Parsing blif format according to the following spec:
// https://course.ece.cmu.edu/~ee760/760docs/blif.pdf
bool BlifParser::parse(std::string& file_contents)
{
  return blif_parser::parse(file_contents.begin(), file_contents.end(), this);
}
}  // namespace rmp
