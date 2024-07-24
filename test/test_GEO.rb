$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'test/unit'
require 'GEO'

class TestGEO < Test::Unit::TestCase
  def _test_query
    assert_equal ["GDS1761"], GEO.query("NCI60[Title]")
  end

  def _test_comparisons
    dataset = "GDS1761"
    info = GEO.dataset_info dataset
    subset_comparisons =  GEO.subset_comparisons info[:subsets]
    assert subset_comparisons['cell line'].include? ['CNS tumor','breast tumor']
  end

  def _test_comparisons_alt
    dataset = "GDS777"
    info = GEO.dataset_info dataset
    subset_comparisons =  GEO.subset_comparisons info[:subsets]
    assert subset_comparisons['cell line'].include? ['CNS tumor','breast tumor']
  end

  def test_codes
    platform = "GPL90"
    sss 0
    GEO[platform].codes.produce
  end
end

