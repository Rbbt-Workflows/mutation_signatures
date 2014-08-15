module Sample

  dep :genomic_mutations
  dep MutationSignatures, :mutation_context, :mutations => :genomic_mutations, :organism => :organism, :watson => :watson
  task :mutation_signature => :tsv do
    TSV.get_stream step(:mutation_context)
  end
end
