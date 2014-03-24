require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/GE'
require 'GEO'

require 'rbbt/matrix'
require 'rbbt/matrix/differential'
require 'rbbt/matrix/barcode'

module GEO
  extend Workflow

  helper :organism do
    "Hsa"
  end

  helper :dataset_info do |dataset|
    GEO[dataset]["info.yaml"].yaml
  end

  helper :identifiers do |dataset|
    if dataset =~ /^GDS|GSE/
      platform = dataset_info(dataset)[:platform]
      GEO[platform].codes
    else
      Organism.identifiers(organism)
    end
  end

  helper :translate do |dataset,keys,format|
    log :translate, "Translating keys to #{ format }"

    if dataset =~ /^GDS|GSE/
      platform = dataset_info(dataset)[:platform]
      codes = GEO[platform].codes
      if TSV.parse_header(codes).all_fields.include? format
        index = codes.index :target => format
      end
    end

    index ||= Organism.identifiers(organism).index :target => format, :persist => true

    keys.collect{|k| index[k] }
  end

  helper :matrix do |dataset|
    if dataset =~ /^GDS|GSE/
      dataset_info = dataset_info dataset

      samples       = dataset_info[:subsets]
      value_type    = dataset_info[:value_type]

      platform = dataset_info[:platform]
      format     = GEO[platform].codes.open{|io| TSV.parse_header(io).key_field } 

      log :geo, "Processing #{ dataset } files from GEO" unless File.exists?(GEO[dataset].values.find)
      data = GEO[dataset].values 

      log :matrix, "Setting matrix for #{ dataset }"
      Matrix.new data.find.produce, samples, value_type, format
    else
      raise "Matrix not found: #{ dataset }" unless File.exists?(dataset)
      format = TSV.parse_header(dataset).key_field

      log :matrix, "Setting matrix for #{ dataset }"
      Matrix.new Path.setup(dataset.dup), nil, nil, format
    end
  end

  helper :matrix_to_gene do |matrix, dataset|
    log :translating, "Translating #{ dataset } matrix to known gene ids"
    begin
      matrix = matrix.to_gene(identifiers(dataset)) 
    rescue FieldNotFoundError
      raise ParameterException, "Could not identify probes in #{ dataset }. Cannot translate to gene -- #{identifiers(dataset).find}"
    end 
  end

  input :query, :string, "Query"
  task :query => :array do |query|
    raise ParameterException, "No query" if query.nil? or query.empty?
    GEO.query query
  end
  export_asynchronous :query
  
  input :dataset, :string, "Dataset code", nil
  task :sample_info => :tsv do |dataset|
    dataset_info = dataset_info dataset

    if dataset_info.include? :sample_info
      sample_info = dataset_info[:sample_info]

      tsv = TSV.setup({}, :key_field => "Sample", :fields => ["Title"], :type => :single)

      sample_info.each do |sample,info|
        tsv[sample] = info[:title]
      end

      tsv
    else
      subsets = dataset_info[:subsets]

      tsv = TSV.setup({}, :key_field => "Sample", :fields => [], :type => :list)
      subsets.each do |factor,values|
        factor_samples = {}
        values.each do |value,samples|
          samples.split(",").each do |sample|
            factor_samples[sample] = value
          end
        end
        factor_samples.each do |sample,value|
          tsv[sample] ||= []
          tsv[sample] << value
        end
        tsv.fields << factor
      end

      tsv
    end
  end
  export_asynchronous :sample_info

  input :dataset, :string, "Dataset code", nil
  input :to_gene, :boolean, "Average probes by gene", false
  task :matrix => :tsv do |dataset,to_gene|
    raise ParameterException, "No dataset provided" if dataset.nil?
    matrix = matrix dataset

    matrix = matrix_to_gene matrix, dataset if to_gene

    FileUtils.cp matrix.data_file, path
  end
  export_asynchronous :matrix

  input :dataset, :string, "Dataset code", nil
  input :main, :string, "Main samples", nil
  input :contrast, :string, "Contrat samples", nil
  input :to_gene, :boolean, "Average probes by gene", false
  task :differential => :tsv do |dataset,main, contrast,to_gene|
    raise ParameterException, "No dataset provided" if dataset.nil?
    matrix = matrix dataset

    matrix = matrix_to_gene matrix, dataset if to_gene

    log :differential, "Running differential for #{ dataset } -- #{Misc.fingerprint [main, contrast]}"
    matrix.differential main, contrast, path
  end
  export_asynchronous :differential

  dep :differential
  input :threshold, :float, "Cuttoff", 0.005
  returns "Ensembl Gene ID"
  task :up_genes => :array do |threshold|
    diff = step(:differential).load
    dataset = step(:differential).info[:inputs][:dataset]
    to_gene = step(:differential).info[:inputs][:to_gene]

    diff.unnamed = true
    keys = diff.select("adjusted.p.values"){|p| p > 0 and p < threshold }.keys

    to_gene ? translate(dataset, keys, "Ensembl Gene ID") : keys
  end
  export_asynchronous :up_genes

  dep :differential
  input :threshold, :float, "Cuttoff", 0.005
  returns "Ensembl Gene ID"
  task :down_genes => :array do |threshold|
    diff = step(:differential).load
    dataset = step(:differential).info[:inputs][:dataset]
    to_gene = step(:differential).info[:inputs][:to_gene]

    diff.unnamed = true
    keys = diff.select("adjusted.p.values"){|p| p < 0 and p.abs < threshold }.keys

    to_gene ? translate(dataset, keys, "Ensembl Gene ID") : keys
  end
  export_asynchronous :down_genes


  input :dataset, :string, "Dataset code", nil
  input :to_gene, :boolean, "Average probes by gene", false
  input :standard_deviations, :float, "Standard deviations from the lower mode", 2
  task :barcode => :tsv do |dataset,to_gene, standard_deviations|
    raise ParameterException, "No dataset provided" if dataset.nil?
    matrix = matrix dataset

    matrix = matrix_to_gene matrix, dataset if to_gene

    matrix.barcode path, standard_deviations
  end
  export_asynchronous :barcode
end
