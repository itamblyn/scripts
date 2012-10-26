#!/usr/bin/perl
#
# gmedtwitpost - set your twitter status in a one line simple perl script with proxy support.
# gmedtwitpost is a simple perl script which allows status updates only - Designed for automated systems to call.
#
# Copyright 2008, Liam Gladdy, http://www.gladdymedia.com/.
#
##############################################################
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.  
#
##############################################################
#
# configuration: uncomment $httpproxy if required, populate $username and $password.
# usage: 'twitter.pl NEW_STATUS'
# version: 0.2, addition of 'source' parameter, as requested by twitter devs
# requires: LWP
#
#

my $username;
my $password;
my $httpproxy = "none";
my $version;

# Uncomment the following line and enter your httpproxy if required, eg: http://squid.proxy.com:1234
#$httpproxy = "http://your.proxyserver.com:port";

# Username and Password
$username = "itamblyn";
$password = "iwanttobelieve";


$version = "0.2";

use strict;
use LWP;

my $useragent = LWP::UserAgent->new;
$useragent->credentials('twitter.com:80','Twitter API',$username,$password);
$useragent->default_headers->push_header("X-Twitter-Client" => "internet");
$useragent->default_headers->push_header("X-Twitter-Client-URL:" => "http://chaffey.phys.dal.ca/~itamblyn/twitter.xml");
$useragent->default_headers->push_header("X-Twitter-Client-Version" => "$version");

if ($httpproxy ne "none") {
$useragent->proxy('http', $httpproxy);
}

if ($#ARGV < 0) {
	print "[FATAL] No Status Given.. What am i supposed to post?\n";
	exit 1;
}

my $doupdate = $useragent->post('http://twitter.com/statuses/update.json', [ "status" => join(' ',@ARGV), "source" => "internet" ]);
die "Invalid response. Twitter returned: ",$doupdate->content."\n" unless $doupdate->content_type =~ m/^application\/json/;
