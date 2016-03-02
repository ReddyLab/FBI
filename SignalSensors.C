/****************************************************************
 SignalSensors.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SignalSensors.H"
using namespace std;
using namespace BOOM;

SignalSensors::SignalSensors()
  : startCodonSensor(NULL), stopCodonSensor(NULL), shortStartSensor(NULL),
    donorSensor(NULL), acceptorSensor(NULL)
{
  //ctor
}



void SignalSensors::setConsensuses()
{
  /*
  SignalSensor *startCodonSensor;
  SignalSensor *shortStartSensor; // just the ATG
  SignalSensor *stopCodonSensor;
  SignalSensor *donorSensor;
  SignalSensor *acceptorSensor;
  Set<String> donorConsensuses;
  Set<String> acceptorConsensuses;
  Set<String> startCodons;
  Set<String> stopCodons;
   */
}



void SignalSensors::setConsensuses(SignalSensor *sensor,
				   const Set<String> &consensuses)
{
  StringMap<char> &known=sensor->getConsensuses();
  for(Set<String>::const_iterator cur=consensuses.begin(), end=
	consensuses.end() ; cur!=end ; ++cur) {
    const String &s=*cur;
    if(!known.isDefined(s)) sensor->addConsensus(s);
  }
}



