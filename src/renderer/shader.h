#pragma once

#include <string>
#include <cstdint>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>

namespace botany {

class Shader {
public:
    Shader() = default;
    bool load(const std::string& vert_path, const std::string& frag_path);
    void use() const;
    void set_mat4(const std::string& name, const glm::mat4& mat) const;
    void set_vec3(const std::string& name, const glm::vec3& vec) const;
    uint32_t id() const { return program_; }

private:
    uint32_t program_ = 0;
};

} // namespace botany
