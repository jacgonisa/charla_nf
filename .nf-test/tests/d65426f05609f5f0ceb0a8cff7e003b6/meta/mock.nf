import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies

include { generate_readmers_kmc  } from '/home/jg2070/Desktop/PhD/crossover/charla_nf/modules/local/tests/../generate_readmers_kmc.nf'


// include test process
include { generate_histogram_kmc } from '/home/jg2070/Desktop/PhD/crossover/charla_nf/modules/local/tests/../generate_histogram_kmc.nf'

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
                    'test_histo',
                    21
                )
                
        generate_readmers_kmc(*input)
    }
    

    // process mapping
    def input = []
    
                input[0] = generate_readmers_kmc.out
                
    //----

    //run process
    generate_histogram_kmc(*input)

    if (generate_histogram_kmc.output){

        // consumes all named output channels and stores items in a json file
        for (def name in generate_histogram_kmc.out.getNames()) {
            serializeChannel(name, generate_histogram_kmc.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = generate_histogram_kmc.out as Object[]
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
