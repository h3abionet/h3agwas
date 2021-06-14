FROM openjdk:latest
COPY . /usr/java
RUN mv /usr/java/beagle /usr/bin/beagle
RUN mv /usr/java/bref3 /usr/bin/bref3
RUN chmod +x /usr/bin/beagle
RUN chmod +x /usr/bin/bref3
WORKDIR /home/data
