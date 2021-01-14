/*
Author: Kenta
Description: Reads XML meta data from tif image and sets channel name
 */

// Import
import groovy.xml.XmlParser

// Destination of XML file
def filename = "/Users/yokote.k/Desktop/test/test.xml"

// Parse in xml
def parser = new XmlParser()
def metadata = parser.parse(filename)

// Extract channel names
def names = metadata.image.channels.channel.@name as String[]

// Remove trailing whitespace
names = names.collect { it -> it.substring(0, it.length()-1)}

// Use GroovyShell.evaluate since setChannelNames does not take String[] as arg
names = names.collect {it -> "'" + it + "'"}
def stringNames = String.join(", ", names)
def command = "import static qupath.lib.gui.scripting.QPEx.*\n"+
                "setChannelNames("+stringNames+")"

// Execute command
new GroovyShell().evaluate(command)

println("Renamed Channels")