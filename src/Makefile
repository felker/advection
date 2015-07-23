CC = gcc
CCFLAGS = 
LINK = gcc
LINKFLAGS = 

all:	advect

%.o:	%.c
	$(CC) $(CCFLAGS) -c $< -o $@

advect:	main.o visit_writer.o	
	$(LINK) $(LINKFLAGS) main.o visit_writer.o -o advect

clean:
	rm *.o advect