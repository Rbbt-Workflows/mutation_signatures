if Module === Sample and Workflow === Sample
  module Sample

    #dep :genomic_mutations
    #dep :organism
    #dep :watson
    #dep MutationSignatures, :context_change_count, :mutations => :genomic_mutations, :organism => :organism, :watson => :watson
    #task :mutation_signature => :tsv do
    #  TSV.get_stream step(:context_change_count)
    #end
    #
    dep :genomic_mutations
    dep :organism
    dep :watson
    dep MutationSignatures, :context_change_count, :mutations => :genomic_mutations, :organism => :organism, :watson => :watson
    task :context_change_count => :tsv do
      TSV.get_stream step(:context_change_count)
    end

    dep :context_change_count
    dep MutationSignatures, :assign_signatures_from_changes, :changes => :context_change_count
    task :mutation_signature => :tsv do
      TSV.get_stream step(:assign_signatures_from_changes)
    end
  end
end
