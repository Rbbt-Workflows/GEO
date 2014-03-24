require 'rbbt'

class Matrix

  class << self
    attr_accessor :matrix_dir
    def matrix_dir
      @matrix_dir ||= Rbbt.var.matrices
    end
  end

  attr_accessor :data_file, :labels, :value_type, :format
  def initialize(data_file, labels, value_type, format)
    @data_file = data_file
    @labels = labels
    @value_type = value_type
    @format = format
  end

  def all_samples
    @all_samples ||= TSV.parse_header(@data_file).fields
  end

  def comparison(subsets, main, contrast)
    if main.index "="
      main_factor, main_value = main.split "=" 
      raise ParameterException, "Main selection not understood" if subsets[main_factor].nil? or subsets[main_factor][main_value].nil?
      main_samples = subsets[main_factor][main_value].split ','
    else
      main_samples = main.split(/[|,\n]/)
    end

    if contrast
      if contrast.index "="
        contrast_factor, contrast_value = contrast.split "=" 
        raise ParameterException, "Contrast selection not understood" if subsets[contrast_factor].nil? or subsets[contrast_factor][contrast_value].nil?
        contrast_samples = subsets[contrast_factor][contrast_value].split ','
      else
        contrast_samples = contrast.split(/[|,\n]/)
      end
    else
      if subsets and defined? main_factor
        contrast_samples = subsets[main_factor].values.collect{|s| s.split ',' }.flatten.uniq - main_samples
      else
        contrast_samples = all_samples - main_samples
      end
    end

    [main_samples, contrast_samples]
  end


  def to_gene(identifiers = nil)
    require 'rbbt/tsv/change_id'
    file = Persist.persist(data_file, :tsv, :prefix => "GENE", :dir => Matrix.matrix_dir, :no_load => true) do
      identifiers = [Organism.identifiers("Hsa"), identifiers].compact

      data_file.tsv(:cast => :to_f).change_key("Ensembl Gene ID", :identifiers => identifiers) do |v|
        Misc.mean(v.compact)
      end
    end
    Matrix.new file, labels, value_type, format
  end
end
