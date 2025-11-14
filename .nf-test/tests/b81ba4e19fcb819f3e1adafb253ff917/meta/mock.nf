import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies

include { generate_readmers_kmc  } from '/home/jg2070/Desktop/PhD/crossover/charla_nf/modules/local/tests/../generate_readmers_kmc.nf'


// include test process
include { get_counts_kmc } from '/home/jg2070/Desktop/PhD/crossover/charla_nf/modules/local/tests/../get_counts_kmc.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()


workflow {

    // run dependencies
    
    {
        def input = []
        
                input[0] = tuple(
                    file('tests/data/test_sample.fasta', checkIfExists: true),
                    'test_sample',
                    21
                )
                
        generate_readmers_kmc(*input)
    }
    

    // process mapping
    def input = []
    
                def kmc_pre = generate_readmers_kmc.out[0][0]
                def kmc_suf = generate_readmers_kmc.out[0][1]
                def fasta = file('tests/data/test_sample.fasta', checkIfExists: true)
                def sample_id = 'test_sample'
                def kmer_size = 21

                input[0] = tuple(kmc_pre, kmc_suf, fasta, sample_id, kmer_size)
                
    //----

    //run process
    get_counts_kmc(*input)

    if (get_counts_kmc.output){

        // consumes all named output channels and stores items in a json file
        for (def name in get_counts_kmc.out.getNames()) {
            serializeChannel(name, get_counts_kmc.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = get_counts_kmc.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
  
}

def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
