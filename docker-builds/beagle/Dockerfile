FROM openjdk:latest
COPY . /usr/java
RUN mv /usr/java/beagle /usr/bin/beagle
RUN mv /usr/java/bref3 /usr/bin/bref3
RUN mv /usr/java/unbref3 /usr/bin/unbref3
RUN mv /usr/java/conform-gt /usr/bin/conform-gt
RUN chmod +x /usr/bin/beagle
RUN chmod +x /usr/bin/bref3
RUN chmod +x /usr/bin/unbref3
RUN chmod +x /usr/bin/conform-gt
WORKDIR /home/data
