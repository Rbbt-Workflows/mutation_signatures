if Module === Sample and Workflow === Sample
  module Sample

    #dep :genomic_mutations
    #dep :organism
    #dep :watson
    #dep MutationSignatures, :context_change_count, :mutations => :genomic_mutations, :organism => :organism, :watson => :watson
    #task :mutation_signature => :tsv do
    #  TSV.get_stream step(:context_change_count)
    #end

    dep :genomic_mutations
    dep :organism
    dep :watson
    dep MutationSignatures, :assign_signatures, :mutations => :genomic_mutations, :organism => :organism, :watson => :watson
    task :mutation_signature => :tsv do
      TSV.get_stream step(:assign_signatures)
    end
  end
end
