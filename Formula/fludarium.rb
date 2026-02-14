class Fludarium < Formula
  desc "Real-time fluid dynamics visualizer"
  homepage "https://github.com/daktu32/fludarium"
  url "https://github.com/daktu32/fludarium/archive/refs/tags/v0.1.0.tar.gz"
  sha256 "e886895dbd4657f203dfa0d5c4f949ee5afbf72880e67c28ca8a1c1edc095605"
  license "MIT"

  depends_on "rust" => :build

  def install
    system "cargo", "install", "--locked", "--root", prefix, "--path", "."
  end

  test do
    assert_match "fludarium", (bin/"fludarium").to_s
  end
end
