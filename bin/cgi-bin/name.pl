#!/usr/bin/perl -w
use strict;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);

print header;
print start_html('Hello');

ReadParse;
print "Here is the form data:<ul>";
    
    foreach $key (keys %in) {
    	print "<li>$key: $in{$key}";
    }
    print "</ul>";

print end_html;
