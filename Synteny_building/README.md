# run.Makeblocks

c++ porting of InferCars makeblocks

- usage 
./run.Makeblocks  config.file >& log

Results and other files are generated in the path where you execute this program

- config file
net/chain directory should be in specific format
\reference species name 
	\target species name
		\chain
			\chr[].chain #chain files split per chromosome
		\net
			\chr[].net #net files split per chromosome

- resolution
resolution for making synteny, minimum synteny block length

- species 
species_name species_tag
 #Species name should be match with the name of chainNet directory
 #Species tag=0: reference species (for only single species)
 #Species tag=1: in-group species (at least 1)
 #Species tag=2: out-group species



