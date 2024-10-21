#! /bin/sh

# run estimation of CO2 emissions of full pipeline fetching job emissions produced for a complete run
curl -X 'GET' \
  'http://tsc-mngd-061.ebi.ac.uk:5000/user/6129202e84824b19ba83d4258684e0b2/footprint/?start=202309142150&stop=202309180950' \
  -H 'accept: application/json'

  
  # Days of GPU runs
  202410160900
  202410161800